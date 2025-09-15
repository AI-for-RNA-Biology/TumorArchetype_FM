# -*- coding: utf-8 -*-
import logging

# Create a logger
logger = logging.getLogger("__main__")
logger.setLevel(logging.DEBUG)  # Set the logging level

# Create a console handler and set the level to debug
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)

# Create a formatter and set it for the handler
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)

# Add the console handler to the logger
logger.addHandler(console_handler)

# Optionally, remove the NullHandler if it's no longer needed
logger.handlers = [h for h in logger.handlers if not isinstance(h, logging.NullHandler)]
