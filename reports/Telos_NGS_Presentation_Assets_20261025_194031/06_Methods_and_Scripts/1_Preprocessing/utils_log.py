import os
import sys
import time
import socket
import getpass
import platform
import logging
from datetime import datetime

def setup_logger(step_name: str, log_dir: str):
    os.makedirs(log_dir, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = os.path.join(log_dir, f"{step_name}_{ts}.log")

    logger = logging.getLogger(step_name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")

    # file handler
    fh = logging.FileHandler(log_path)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    # stdout handler
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    # header
    logger.info("============================================================")
    logger.info(f"STEP: {step_name}")
    logger.info(f"LOG:  {log_path}")
    logger.info(f"USER: {getpass.getuser()}  HOST: {socket.gethostname()}")
    logger.info(f"SYS:  {platform.platform()}")
    logger.info("============================================================")

    return logger, log_path

def write_done_marker(done_path: str, content: str):
    os.makedirs(os.path.dirname(done_path), exist_ok=True)
    with open(done_path, "w") as f:
        f.write(content.rstrip() + "\n")