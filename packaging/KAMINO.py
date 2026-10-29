"""The program PyInstaller turns into the downloadable app (see KAMINO.spec)."""
import multiprocessing

if __name__ == '__main__':
    multiprocessing.freeze_support()     # the measuring worker processes start through this same program
    from kneemeasure.desktop import main
    main()
