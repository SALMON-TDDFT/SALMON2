import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from salmonpy.errors import InputError, RuntimeFailure
from salmonpy.pseudopotentials import PseudopotentialReference, download, silicon_reference, _HTTPSRedirect


class Response(io.BytesIO):
    def geturl(self):
        return 'https://example.org/Si_rps.dat'


class PseudopotentialTests(unittest.TestCase):
    def reference(self, content=b'potential'):
        return PseudopotentialReference('Si_rps.dat', 14, 'https://example.org/Si_rps.dat',
                                       hashlib.sha256(content).hexdigest(), 'Apache-2.0', 'pinned upstream')

    def test_download_permission_identity_and_existing_destination(self):
        with tempfile.TemporaryDirectory() as tmp, patch('salmonpy.pseudopotentials.build_opener') as opener:
            target = Path(tmp)/'new'
            with self.assertRaises(InputError):
                download(self.reference(), target)
            opener.assert_not_called()
            opener.return_value.open.return_value = Response(b'potential')
            path = download(self.reference(), target, allow_download=True)
            self.assertEqual(path.read_bytes(), b'potential')
            receipt = json.loads((target/'download.json').read_text())
            self.assertEqual(receipt['sha256'], self.reference().sha256)
            with self.assertRaises(FileExistsError):
                download(self.reference(), target, allow_download=True)
            self.assertEqual(opener.call_count, 1)

    def test_rejected_download_does_not_publish_directory(self):
        for content, url in ((b'wrong', 'https://example.org/Si_rps.dat'),
                             (b'potential', 'http://example.org/Si_rps.dat')):
            with self.subTest(url=url, content=content), tempfile.TemporaryDirectory() as tmp, patch('salmonpy.pseudopotentials.build_opener') as opener:
                response = Response(content)
                response.geturl = lambda: url
                opener.return_value.open.return_value = response
                target = Path(tmp)/'new'
                with self.assertRaises((RuntimeFailure, InputError)):
                    download(self.reference(), target, allow_download=True)
                self.assertFalse(target.exists())
                self.assertEqual(list(Path(tmp).iterdir()), [])
        with patch('salmonpy.pseudopotentials._MAX_BYTES', 3), tempfile.TemporaryDirectory() as tmp, patch('salmonpy.pseudopotentials.build_opener') as opener:
            opener.return_value.open.return_value = Response(b'potential')
            with self.assertRaisesRegex(RuntimeFailure, 'limit'):
                download(self.reference(), Path(tmp)/'new', allow_download=True)

    def test_reference_validation_and_https_redirect(self):
        fields = vars(self.reference())
        for changes in ({'filename': '../escape'}, {'filename': 'download.json'}, {'filename': 'DOWNLOAD.JSON'}, {'sha256': '0'},
                        {'url': 'http://example.org/p'}, {'url': 'https://user:secret@example.org/p'},
                        {'license': ''}, {'atomic_number': True}):
            with self.subTest(changes=changes), self.assertRaises(InputError):
                PseudopotentialReference(**(fields | changes))
        with self.assertRaises(InputError):
            _HTTPSRedirect().redirect_request(None, None, 302, '', {}, 'http://example.org/p')
        reference = silicon_reference()
        content = Path(__file__).parents[1]/'examples/si_intro_data/Si_rps.dat'
        self.assertEqual(hashlib.sha256(content.read_bytes()).hexdigest(), reference.sha256)
