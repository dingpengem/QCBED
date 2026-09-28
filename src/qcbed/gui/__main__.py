"""Launch the optional configuration editor: python -m qcbed.gui [config.yml]."""
import argparse
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(description='QCBED canonical configuration editor and refinement runner')
    parser.add_argument('config', nargs='?', help='Canonical YAML configuration to open')
    args = parser.parse_args(argv)
    try:
        from PySide6.QtWidgets import QApplication, QMessageBox
        from .window import ConfigWindow
    except ModuleNotFoundError as exc:
        if exc.name and exc.name.startswith('PySide6'):
            parser.exit(2, 'The optional GUI dependency is missing. Install with: pip install -e ".[gui]"\n')
        raise
    app = QApplication(sys.argv[:1])
    window = ConfigWindow()
    if args.config:
        try: window.open_document(args.config)
        except (OSError, ValueError, TypeError) as exc:
            QMessageBox.warning(window, 'Cannot open configuration', str(exc))
    window.show()
    return app.exec()


if __name__ == '__main__':
    raise SystemExit(main())
