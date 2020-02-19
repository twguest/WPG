import logging
import time


def config_log(filename=None):

    date = time.strftime("%x").replace("/", "_")  # date

    if filename == None:
        filename = date

    try:
        logging.basicConfig(
            filename="./out/logs/twg/logs/" + filename + ".log",
            level=logging.DEBUG,
            format="%(asctime)s %(message)s",
            datefmt="%m/%d/%Y %I:%M:%S %p",
        )
    except:
        filepath = r"./out/logs/twg/logs/"
        logging.basicConfig(
            filename=filepath + filename + ".log",
            level=logging.DEBUG,
            format="%(asctime)s %(message)s",
            datefmt="%m/%d/%Y %I:%M:%S %p",
        )

    log = logging.getLogger(__name__)

    return log


def get_log():
    return logging.getLogger(__name__)
