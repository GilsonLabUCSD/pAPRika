import logging
import os

# This is modeled after YANK's nice logging facility.
# https://github.com/choderalab/yank/blob/4dfcc8e127c51c20180fe6caeb49fcb1f21730c6/Yank/utils.py#L78


def config_root_logger(verbose, log_file_path=None):
    """
    Setup the the root logger's configuration.
    The log messages are printed in the terminal and saved in the file specified
    by log_file_path (if not `None`) and printed. Note that logging use sys.stdout
    to print ``logging.INFO`` messages, and stderr for the others. The root logger's
    configuration is inherited by the loggers created by ``logging.getLogger(name)``.
    Different formats are used to display messages on the terminal and on the log
    file. For example, in the log file every entry has a timestamp which does not
    appear in the terminal. Moreover, the log file always shows the module that
    generate the message, while in the terminal this happens only for messages
    of level `WARNING` and higher.

    Parameters
    ----------
    verbose : bool
        Control the verbosity of the messages printed in the terminal. The logger
        displays messages of level ``logging.INFO`` and higher when ``verbose=False``.
        Otherwise those of level ``logging.DEBUG`` and higher are printed.
    log_file_path : str, optional, default = None
        If not `None`, this is the path where all the logger's messages of level
        ``logging.DEBUG`` or higher are saved.
    """

    class TerminalFormatter(logging.Formatter):
        """
        Simplified format for INFO and DEBUG level log messages.
        This allows to keep the logging.info() and debug() format separated from
        the other levels where more information may be needed. For example, for
        warning and error messages it is convenient to know also the module that
        generates them.
        """

        # This is the cleanest way I found to make the code compatible with both
        # Python 2 and Python 3
        simple_fmt = logging.Formatter("%(asctime)-15s: %(message)s")
        default_fmt = logging.Formatter(
            "%(asctime)-15s: %(levelname)s - %(name)s - %(message)s"
        )

        def format(self, record):
            if record.levelno <= logging.INFO:
                return self.simple_fmt.format(record)
            else:
                return self.default_fmt.format(record)

    # Add handler for stdout and stderr messages
    terminal_handler = logging.StreamHandler()
    terminal_handler.setFormatter(TerminalFormatter())
    if verbose:
        terminal_handler.setLevel(logging.DEBUG)
    else:
        terminal_handler.setLevel(logging.INFO)
    logging.root.addHandler(terminal_handler)

    # Add file handler to root logger
    file_format = "%(asctime)s - %(levelname)s - %(name)s - %(message)s"
    if log_file_path is not None:
        file_handler = logging.FileHandler(log_file_path)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(file_format))
        logging.root.addHandler(file_handler)

    # Do not handle logging.DEBUG at all if unnecessary
    if log_file_path is not None:
        logging.root.setLevel(logging.DEBUG)
    else:
        logging.root.setLevel(terminal_handler.level)

    # Setup critical logger file if a logfile is specified
    # No need to worry about MPI due to it already being set above
    if log_file_path is not None:
        basepath, ext = os.path.splitext(log_file_path)
        critical_log_path = basepath + "_CRITICAL" + ext
        # Create the critical file handler to only create the file IF a critical message is sent
        critical_file_handler = logging.FileHandler(critical_log_path, delay=True)
        critical_file_handler.setLevel(logging.CRITICAL)
        # Add blank lines to space out critical errors
        critical_file_format = file_format + "\n\n\n"
        critical_file_handler.setFormatter(logging.Formatter(critical_file_format))
        logging.root.addHandler(critical_file_handler)
