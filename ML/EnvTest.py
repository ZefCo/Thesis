print("Starting Env Test")

try:
    import os
except Exception as e:
    print("Unable to import module os")
    fail_state = True

try:
    import shutil
except Exception as e:
    print("Unable to import module shutil")
    fail_state = True

try:
    import pandas
except Exception as e:
    print("Unable to import module pandas")
    fail_state = True

try:
    import pathlib
except Exception as e:
    print("Unable to import module pathlib")
    fail_state = True

try:
    import re
except Exception as e:
    print("Unable to import module re")
    fail_state = True

try:
    import sys
except Exception as e:
    print("Unable to import module sys")
    fail_state = True

try:
    import tensorflow as tf
    print("Unable to import module tensorflow")
    fail_state = True
    print(tf.__version__)

except Exception as e:
    print("Unable to import module")
    fail_state = True

try:
    import random
except Exception as e:
    print("Unable to import module random")
    fail_state = True

try:
    import datetime
except Exception as e:
    print("Unable to import module datetime")
    fail_state = True

try:
    from matplotlib import pyplot as plt
except Exception as e:
    print("Unable to import module matplotlib")
    fail_state = True

try:
    from PIL import Image
except Exception as e:
    print("Unable to import module PIL")
    fail_state = True


if fail_state:
    print("Evn test FAILED")
else:
    print("Env Test complete: passed")
# if "tf" not in sys.modules:
#     print("TensorFlow not loaded")
# else:
#     print(tf.__version__)
