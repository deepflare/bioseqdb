"""
Logging utilities
@Piotr Styczyński 2021
"""
import inspect
import logging
import traceback
from tqdm import tqdm
from rich.logging import RichHandler

from rich.console import Console

def pretty_print_to_str(*args, **kwargs):
    console = Console()
    with console.capture() as capture:
        console.print(*args, **kwargs)
    str_output = capture.get()
    return str_output[:-1]

def escape_password(password: str) -> str:
    if len(password) < 6:
        return "*" * len(password)
    return password[:1] + ("*" * (len(password) - 3)) + password[-1]


class UniversalLoggerSet:
    """
    Collection of loggers.
    """

    def __init__(self):
        self._loggers = dict()

    def get_logger(self, name: str, logger_type=RichHandler):
        """
        Get universal logger for a given name.
        If the logger does not exist yet, create it and add to the UniversalLoggerSet
        If it exists already, just return it
        :param name: <str> Name of the logger
        :return: Logger instance
        """
        if name in self._loggers:
            return self._loggers[name]

        logger = logging.getLogger(name)
        # RichHandler selected
        if logger_type == RichHandler:
            # initialize RichHandler
            handler = RichHandler(
                level=logging.DEBUG,  # default NOTSET
                show_time=False,  # default true
                show_level=True,
                enable_link_path=True,  # default true
                rich_tracebacks=True,
            )
        # StreamHandler selected
        elif logger_type == logging.StreamHandler:
            handler = logging.StreamHandler()
        logger.addHandler(handler)
        logger.setLevel(logging.DEBUG)

        universal_logger = UniversalLogger(logger, name)
        self._loggers[name] = universal_logger
        return universal_logger

    def __getattr__(self, attr):
        return self.get_logger(attr)


class UniversalLogger:
    """
    Wrapper for Python logging library logger and set of useful logging utilities.
    """

    def __init__(self, logger_instance, name):
        self._logger_instance = logger_instance
        self._name = name

    def _preformat_message(self, message, level):
        """
        Format message before passing to the inner logger.
        :param message: Input message string
        :param level: Level of logging for that message
        :return: Preformatted message
        """
        frame = inspect.currentframe()
        this_frame = frame  # Save current frame.

        frame_no = 0
        while frame.f_back:
            frame = frame.f_back
            if frame_no == 1:
                this_frame = frame
            frame_no = frame_no + 1

        this_frame_info = inspect.getframeinfo(this_frame)
        filename = "/".join(this_frame_info.filename.split("/")[-1:])
        lineno = this_frame_info.lineno
        name_block = "%-10s" % self._name
        source_block = "%-15s |" % f"{filename}:{lineno}"
        return f"{name_block} {source_block} {message}"

    def error(self, message):
        """
        Log error message.
        :param message: Message to be logged
        """
        message: str = f"❌ {message}"
        self._logger_instance.error(self._preformat_message(message, "ERROR"))

    def warn(self, message):
        """
        Log warning message.
        :param message: Message to be logged
        """
        message: str = f"🚸 {message}"
        self._logger_instance.warn(self._preformat_message(message, "WARN"))

    def ok(self, message):
        """
        Log info message.
        :param message: Message to be logged
        """
        message: str = f"✅ {message}"
        self._logger_instance.info(self._preformat_message(message, "INFO"))

    def info(self, message):
        """
        Log info message.
        :param message: Message to be logged
        """
        message: str = f" {message}"
        self._logger_instance.info(self._preformat_message(message, "INFO"))

    def debug(self, message):
        """
        Log debug message.
        :param message: Message to be logged
        """
        message: str = f"❗ {message}"
        self._logger_instance.debug(self._preformat_message(message, "DEBUG"))


def with_logging(logger: UniversalLogger):
    def decorator(function):
        def wrapper(*args, **kwargs):
            try:
                result = function(*args, **kwargs)
                return result
            except Exception as e:
                logger.error(
                    f"Function {function.__name__} raised an exception. Exception details: {e}"
                )
                logger.error(f"Stack trace: {traceback.format_exc()}")

        return wrapper

    return decorator


LOGS = UniversalLoggerSet()

if __name__ == "__main__":
    LOGS.task.info("I am INFO")
    LOGS.task.error("I am ERROR")
    LOGS.task.debug("I am DEBUG")