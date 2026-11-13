"""Build the same CMake LD target used by wheels, into a source checkout."""
from pathlib import Path
import shutil
import subprocess
import tempfile


def main():
    root = Path(__file__).resolve().parent.parent
    with tempfile.TemporaryDirectory(prefix='salmonpy-ld-') as tmp:
        tmp = Path(tmp)
        subprocess.run(['cmake', '-S', str(root / 'packaging/ld'), '-B', str(tmp / 'build'),
                        '-DCMAKE_BUILD_TYPE=Release', '-DCMAKE_INSTALL_PREFIX=' + str(tmp / 'install')], check=True)
        subprocess.run(['cmake', '--build', str(tmp / 'build'), '--config', 'Release'], check=True)
        subprocess.run(['cmake', '--install', str(tmp / 'build'), '--config', 'Release'], check=True)
        target = root / 'salmonpy/_native'
        target.mkdir(exist_ok=True)
        for source in (tmp / 'install/salmonpy/_native').iterdir():
            shutil.copy2(source, target / source.name)
        print(target)

if __name__ == '__main__':
    main()
