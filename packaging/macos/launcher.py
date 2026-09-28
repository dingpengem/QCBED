"""Frozen-only entry: dispatch spawn workers before importing Qt or qcbed."""
if __name__ == '__main__':
    import multiprocessing
    multiprocessing.freeze_support()

    import sys
    if len(sys.argv) == 3 and sys.argv[1] == '--smoke-test':
        from app_smoke import main
        raise SystemExit(main(sys.argv[2]))

    from qcbed.gui.__main__ import main
    raise SystemExit(main())
