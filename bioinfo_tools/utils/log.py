import logging
import sys
import os


class Log:
    def __init__(self, logger_name=None, log_level='INFO', log_dir=None, logger:logging.Logger=None, **kwargs):
        self.logger:logging.Logger = None
        self.logger_name = logger_name
        self.log_level = log_level
        self.log_dir = log_dir
        
        if logger and isinstance(logger, logging.Logger):
            self.logger = logger
        elif self.logger_name:
            # retrieve existing logger if already set
            if self.logger_name in logging.Logger.manager.loggerDict.keys():
                self.logger = logging.getLogger(self.logger_name)
            
            # set filehandler logger
            elif self.log_dir:
                if not os.path.exists(log_dir):
                    os.makedirs(log_dir, exist_ok=True)
                log_filepath = os.path.join(self.log_dir, self.logger_name + '.log')
                self.logger = self.get_file_logger(filepath=log_filepath)
            
            # set stream stdout logger
            else:
                self.logger = self.get_stream_logger()

    def get_file_logger(self, filepath) -> logging.Logger:
        handler_formatter = logging.Formatter(fmt='%(message)s')

        handler = logging.FileHandler(filename=filepath, mode="w", encoding='utf8')
        handler.setFormatter(handler_formatter)
        handler.setLevel(self.log_level)

        logger = logging.getLogger(self.logger_name)
        logger.setLevel(self.log_level)
        logger.addHandler(handler)

        return logger

    def get_stream_logger(self, fmt='%(asctime)s  %(message)s') -> logging.Logger:
        handler_formatter = logging.Formatter(fmt=fmt, datefmt='[%Y/%m/%d %I:%M:%S]')

        handler = logging.StreamHandler(stream=sys.stdout)
        handler.setFormatter(handler_formatter)
        handler.setLevel(self.log_level)

        logger = logging.getLogger(self.logger_name)
        logger.setLevel(self.log_level)
        logger.addHandler(handler)

        return logger

    def log(self, *msg, **kwargs):
        level = kwargs.get('level', 'INFO')
        
        # convert level to code if nedeed
        if isinstance(level, str):
            level = getattr(logging, level)

        if self.logger:
            self.logger.log(level, " ".join(msg))
        else:
            print(*msg)
