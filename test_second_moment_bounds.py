"""Exact tests of returned atoms in the ORIGINAL longitudinal/transverse D.

Quadratic complex arithmetic uses t^2=-gamma_squared, so even irrational
growth rates have exact rational real and t coefficients. No floating point
root finding, denominator clearing, or reduced scalar equation is used.
"""

from fractions import Fraction as F
import unittest

from second_moment_bounds import independent_r_bounds, second_moment_bounds


class QC:
    """a+b*t in Q[t]/(t^2+g2), with g2 strictly positive."""
    def __init__(self, a, b=0, g2=1):
        self.a, self.b, self.g2 = map(F, (a, b, g2))

    def convert(self, value):
        if isinstance(value, QC):
            if value.g2 != self.g2:
                raise ValueError("quadratic fields must have the same gamma_squared")
            return value
        return QC(value, 0, self.g2)

    def __add__(self, value):
        v = self.convert(value)
        return QC(self.a + v.a, self.b + v.b, self.g2)

    __radd__ = __add__

    def __neg__(self):
        return QC(-self.a, -self.b, self.g2)

    def __sub__(self, value):
        return self + -self.convert(value)

    def __rsub__(self, value):
        return self.convert(value) + -self

    def __mul__(self, value):
        v = self.convert(value)
        return QC(self.a*v.a-self.b*v.b*self.g2,
                  self.a*v.b+self.b*v.a, self.g2)

    __rmul__ = __mul__

    def inverse(self):
        d = self.a*self.a + self.b*self.b*self.g2
        if d == 0:
            raise ZeroDivisionError
        return QC(self.a/d, -self.b/d, self.g2)

    def __truediv__(self, value):
        return self * self.convert(value).inverse()

    def __rtruediv__(self, value):
        return self.convert(value) * self.inverse()


def original_dispersion(attainer):
    w = QC(attainer.omega_real, 1, attainer.gamma_squared)
    integrals = [QC(0, 0, attainer.gamma_squared) for _ in range(3)]
    for sign, atoms in ((1, attainer.plus_atoms), (-1, attainer.minus_atoms)):
        for atom in atoms:
            denominator = w - attainer.k * atom.u
            for j in range(3):
                integrals[j] += sign * atom.weight * atom.u**j / denominator
    i0, i1, i2 = integrals
    d = (1+i0)*(-1+i2)-i1*i1
    t = -1+(i0-i2)/2
    return d, t


class SecondMomentTest(unittest.TestCase):
    def test_exact_original_dispersion_and_positive_moments(self):
        # 225 rational moment pairs at three number scales, both branches.
        nodes = sorted(set(F(j, 8) for j in range(9)) |
                       {F(1, 100), F(1, 3), F(2, 3), F(99, 100), F(1, 7), F(6, 7)})
        cases = 0
        for n in (F(1, 3), F(1), F(7, 2)):
            for ap in nodes:
                for am in nodes:
                    result = second_moment_bounds(n, n*ap, n*am)
                    if result.stable_reason:
                        continue
                    for name in ("longitudinal", "transverse"):
                        branch = getattr(result, name)
                        att = branch.attainer
                        self.assertGreater(att.k, 0)
                        self.assertGreater(branch.gamma_squared, 0)
                        self.assertEqual(att.gamma_squared, branch.gamma_squared)
                        for atoms, q in ((att.plus_atoms, n*ap), (att.minus_atoms, n*am)):
                            self.assertEqual(sum(a.weight for a in atoms), n)
                            self.assertEqual(sum(a.weight*a.u for a in atoms), 0)
                            self.assertEqual(sum(a.weight*a.u*a.u for a in atoms), q)
                            self.assertTrue(all(a.weight > 0 and -1 <= a.u <= 1 for a in atoms))
                            self.assertEqual(len(set(a.u for a in atoms)), len(atoms))
                        d, t = original_dispersion(att)
                        residual = d if name == "longitudinal" else t
                        self.assertEqual((residual.a, residual.b), (0, 0))
                        cases += 1
        self.assertEqual(cases, 1170)

    def test_stable_controls(self):
        for qp, qm in ((0, 0), (1, 0), (0, 1), (1, 1), (1, F(1, 3)), (F(2, 3), 1)):
            result = second_moment_bounds(1, qp, qm)
            self.assertIsNotNone(result.stable_reason)
            self.assertEqual(result.global_gamma_squared, 0)
            for branch in (result.longitudinal, result.transverse):
                self.assertEqual(branch.gamma_squared, 0)
                self.assertIsNone(branch.attainer)

    def test_zero_central_node(self):
        result = second_moment_bounds(1, 0, F(1, 2))
        for branch in (result.longitudinal, result.transverse):
            self.assertEqual(len(branch.attainer.plus_atoms), 1)
            atom = branch.attainer.plus_atoms[0]
            self.assertEqual((atom.u, atom.weight), (0, 1))

    def test_r_interface_and_fraction_strings(self):
        self.assertEqual(independent_r_bounds("3/2", "1/3", "5/6"),
                         second_moment_bounds("3/2", "7/6", "2/3"))
        self.assertEqual(second_moment_bounds("1", "0.2", "0.4").global_gamma_squared,
                         F(12, 25))

    def test_invalid_inputs(self):
        for values in ((0, 0, 0), (-1, 0, 0), (1, -1, 0), (1, 0, 2)):
            with self.assertRaises(ValueError):
                second_moment_bounds(*values)
        with self.assertRaises(ValueError):
            independent_r_bounds(1, 2, 0)
        with self.assertRaises(ValueError):
            QC(0, g2=1) + QC(0, g2=2)


if __name__ == "__main__":
    unittest.main()
