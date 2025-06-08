# SNPsnap.spec

block_cipher = None

# Add your icon file to datas
a = Analysis(
    ['pyqt.py'],
    pathex=['.'],
    binaries=[],
    datas=[
        ('mafft.bat', '.'),                       # Copies project_root/mafft.bat
        ('usr', 'usr'),                           # Copies project_root/usr
        ('app.ico', '.')                          
    ],
    hiddenimports=[
        'xlsxwriter',
        'Bio.AlignIO',
        'Bio.SeqIO',
        'Bio.SeqRecord',
        'Bio.Seq',
        'Bio.Data',
        'pandas',
        'numpy',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(
    a.pure,
    a.zipped_data,
    cipher=block_cipher
)

exe_options = {
    'name': 'SNPsnap',
    'debug': False,
    'bootloader_ignore_signals': False,
    'strip': False,
    'upx': True,
    'upx_exclude': [],
    'runtime_tmpdir': None,
    'console': False,                
    'disable_windowed_traceback': False,
    'argv_emulation': False,
    'target_arch': None,
    'codesign_identity': None,
    'entitlements_file': None,
    'icon': 'app.ico',               
}

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    **exe_options
)

# Build the one-file executable directly from this spec
# Run:
#   pyinstaller SNPsnap.spec
