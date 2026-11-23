"""Explicit, hash-verified pseudopotential acquisition; never used implicitly."""
import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import tempfile
from urllib.parse import urlsplit
from urllib.request import HTTPRedirectHandler, build_opener

from .errors import InputError, RuntimeFailure

_MAX_BYTES = 64 * 1024 * 1024


def _https(url):
    parts = urlsplit(url)
    if parts.scheme != 'https' or not parts.hostname or parts.username is not None or parts.password is not None or parts.fragment:
        raise InputError('pseudopotential URL must be HTTPS without credentials or fragments')
    return url


class _HTTPSRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        _https(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


@dataclass(frozen=True)
class PseudopotentialReference:
    """Pinned file identity and declared rights, not a physical validation claim."""
    filename: str
    atomic_number: int
    url: str
    sha256: str
    license: str
    source: str
    catalog_id: str | None = None
    provider: str | None = None
    version: str | None = None
    element: str | None = None
    xc: str | None = None
    relativity: str | None = None
    valence_electrons: int | None = None
    nlcc: bool | None = None
    has_so: bool | None = None
    accuracy: str | None = None
    upstream_commit: str | None = None
    file_format: str | None = None
    semicore_shells: tuple | None = None
    valence_configuration: str | None = None

    def __post_init__(self):
        if (not self.filename or self.filename.casefold() in ('.', '..', 'download.json') or
                Path(self.filename).name != self.filename or
                any(c in self.filename for c in '/\\:') or
                any(ord(c) < 32 for c in self.filename)):
            raise InputError('filename must be one ordinary file name')
        if type(self.atomic_number) is not int or not 1 <= self.atomic_number <= 118:
            raise InputError('atomic_number must be an integer from 1 to 118')
        if self.semicore_shells is not None:
            if not isinstance(self.semicore_shells,(list,tuple)) or any(not isinstance(x,str) or not x for x in self.semicore_shells):
                raise InputError('semicore_shells must contain shell names, or be unknown')
            object.__setattr__(self,'semicore_shells',tuple(self.semicore_shells))
        if self.relativity not in (None, 'scalar-relativistic', 'fully-relativistic'):
            raise InputError('unknown relativity; no implicit scalar/SOC conversion')
        for flag in (self.has_so, self.nlcc):
            if flag is not None and type(flag) is not bool:
                raise InputError('has_so/nlcc must be boolean or unknown')
        if self.valence_electrons is not None and (type(self.valence_electrons) is not int or
                not 1 <= self.valence_electrons <= self.atomic_number):
            raise InputError('invalid valence_electrons')
        if self.relativity is not None and self.has_so is not None and self.has_so != (self.relativity == 'fully-relativistic'):
            raise InputError('relativity and has_so disagree')
        _https(self.url)
        if len(self.sha256) != 64 or any(c not in '0123456789abcdef' for c in self.sha256):
            raise InputError('sha256 must contain 64 lowercase hexadecimal characters')
        if not self.license.strip() or not self.source.strip():
            raise InputError('license and source are required; verify permission before downloading')


def list_references(*, element=None, xc=None, relativity=None, provider=None):
    """Search only the installed catalog. No network, substitution or conversion."""
    data = json.loads(Path(__file__).with_name('pseudopotential_catalog.json').read_text(encoding='utf-8'))
    if data.get('schema') != 'salmonpy.pseudopotential-catalog/1':
        raise InputError('unsupported pseudopotential catalog schema')
    selected = []
    for item in data['references']:
        reference = PseudopotentialReference(**item)
        conditions = ((element, reference.element), (xc, reference.xc),
                      (relativity, reference.relativity), (provider, reference.provider))
        if all(wanted is None or (actual is not None and wanted.casefold() == actual.casefold())
               for wanted, actual in conditions):
            selected.append(reference)
    return tuple(selected)


def catalog_reference(identifier):
    """Select an exact catalog ID; unavailable choices never fall back."""
    matches = [r for r in list_references() if r.catalog_id == identifier]
    if len(matches) != 1:
        raise InputError('unknown or nonunique pseudopotential catalog ID: '+str(identifier))
    return matches[0]


def _check_upf(reference, data):
    import xml.etree.ElementTree as ET
    try:
        header = ET.fromstring(data).find('PP_HEADER')
        if header is None or header.get('pseudo_type', '').upper() != 'NC':
            raise InputError('catalog requires an NC UPF header')
        values = header.attrib
        flags = lambda key: values.get(key, 'F').strip().lower() in ('t','true','.true.')
        if flags('is_paw') or flags('is_ultrasoft'):
            raise InputError('catalog requires norm-conserving data')
        expected = [(reference.element, values.get('element', '').strip()),
                    (reference.xc, values.get('functional', '').strip()),
                    (reference.has_so, flags('has_so')),
                    (reference.nlcc, flags('core_correction'))]
        if reference.valence_electrons is not None:
            expected.append((reference.valence_electrons, float(values['z_valence'])))
        if reference.semicore_shells is not None:
            shells={e.get('label','').lower().strip() for e in ET.fromstring(data).findall('PP_PSWFC/*')}
            if not set(reference.semicore_shells)<=shells:
                raise InputError('UPF lacks a declared semicore shell')
        if any(wanted is not None and wanted != actual for wanted,actual in expected):
            raise InputError('UPF header differs from pinned catalog metadata')
    except (ET.ParseError, KeyError, ValueError) as exc:
        raise InputError('invalid catalog UPF metadata') from exc


def silicon_reference():
    """Return the pinned Apache-2.0 Si example reference, without network access."""
    data = Path(__file__).with_name('examples') / 'si_intro_data' / 'provenance.json'
    provenance = json.loads(data.read_text(encoding='utf-8'))
    item = provenance['pseudopotential']
    commit = provenance['upstream_commit']
    return PseudopotentialReference(
        'Si_rps.dat', 14,
        f"https://raw.githubusercontent.com/SALMON-TDDFT/SALMON2/{commit}/{item['path']}",
        item['sha256'], item['license'], item['source_url'])


def download(reference, directory, *, allow_download=False):
    """Download one authorized reference into a NEW directory and record provenance.

    The caller must establish permission under the declared license. No parsing,
    conversion, compiler invocation or calculation is performed. A changed hash
    is an error; it never updates the reference or replaces an existing file.
    """
    if not isinstance(reference, PseudopotentialReference):
        raise InputError('reference must be a PseudopotentialReference')
    if allow_download is not True:
        raise InputError('explicit allow_download=True is required')
    destination = Path(directory).expanduser().absolute()
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    if not destination.parent.is_dir():
        raise InputError('destination parent directory must already exist')
    with tempfile.TemporaryDirectory(prefix='.salmonpy-pseudo-', dir=destination.parent) as tmp:
        staged = Path(tmp) / reference.filename
        digest = hashlib.sha256()
        size = 0
        with build_opener(_HTTPSRedirect()).open(reference.url, timeout=60) as response, staged.open('xb') as out:
            _https(response.geturl())
            while True:
                block = response.read(min(1024 * 1024, _MAX_BYTES + 1 - size))
                if not block:
                    break
                size += len(block)
                if size > _MAX_BYTES:
                    raise RuntimeFailure('pseudopotential exceeds 64 MiB download limit')
                digest.update(block)
                out.write(block)
        if not size or digest.hexdigest() != reference.sha256:
            raise RuntimeFailure('pseudopotential checksum mismatch or empty file')
        if reference.file_format == 'UPF2':
            _check_upf(reference, staged.read_bytes())
        receipt = {'schema': 'salmonpy.pseudopotential-download/1',
                   'reference': asdict(reference), 'bytes': size,
                   'sha256': digest.hexdigest(), 'final_url': response.geturl(),
                   'scope': 'file identity and declared provenance; not physical validation'}
        # Reserve the new directory only after verification. Exclusive creation
        # also rejects a racing creator or symlink without touching their files.
        destination.mkdir()
        with (destination / reference.filename).open('xb') as out:
            out.write(staged.read_bytes())
        with (destination / 'download.json').open('x', encoding='utf-8') as out:
            out.write(json.dumps(receipt, indent=2)+'\n')
    return destination / reference.filename


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', nargs='?', help='new output directory; parent must exist')
    refs = parser.add_mutually_exclusive_group()
    refs.add_argument('--silicon-example', action='store_true')
    refs.add_argument('--reference', help='local JSON containing PseudopotentialReference fields')
    refs.add_argument('--catalog-id', help='exact ID printed by --list')
    parser.add_argument('--list', action='store_true', help='search local metadata only; no download')
    parser.add_argument('--element')
    parser.add_argument('--xc')
    parser.add_argument('--relativity', choices=('scalar-relativistic','fully-relativistic'))
    parser.add_argument('--allow-download', action='store_true', help='explicit authorization to fetch the selected file')
    args = parser.parse_args()
    if args.list:
        if args.directory or args.catalog_id or args.reference or args.silicon_example or args.allow_download:
            parser.error('--list is metadata-only and cannot download')
        print(json.dumps([asdict(r) for r in list_references(element=args.element,xc=args.xc,
                         relativity=args.relativity)],indent=2))
        return
    if not args.directory or not (args.catalog_id or args.reference or args.silicon_example):
        parser.error('download requires a new directory and an exact reference')
    if args.element or args.xc or args.relativity:
        parser.error('filters apply only to --list; download uses an exact reference')
    if args.catalog_id:
        reference = catalog_reference(args.catalog_id)
    elif args.silicon_example:
        reference = silicon_reference()
    else:
        reference = PseudopotentialReference(**json.loads(Path(args.reference).read_text(encoding='utf-8')))
    print(download(reference, args.directory, allow_download=args.allow_download))


if __name__ == '__main__':
    main()
