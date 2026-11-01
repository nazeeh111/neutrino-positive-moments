"""Verify source/data bytes against the delivered manifest; standard library."""
from pathlib import Path
import hashlib,json

def main():
    root=Path(__file__).resolve().parent
    manifest=json.loads((root/'BUNDLE_MANIFEST.json').read_text())
    for name,record in manifest['files'].items():
        relative=Path(name)
        if relative.is_absolute() or '..' in relative.parts:raise ValueError('Unsafe manifest path')
        target=root/relative
        data=target.read_bytes()
        if len(data)!=record['bytes'] or hashlib.sha256(data).hexdigest()!=record['sha256']:
            raise ValueError(f'Manifest mismatch: {name}')
    print(f"Verified {len(manifest['files'])} file hashes and sizes")
if __name__=='__main__':main()
