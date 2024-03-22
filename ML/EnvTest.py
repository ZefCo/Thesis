print("Starting Env Test")

try:
    import os
except Exception as e:
    print("Unable to import os")
    fail_state = True

try:
    import shutil
except Exception as e:
    print("Unable to import shutil")
    fail_state = True

try:
    import pandas
except Exception as e:
    print("Unable to import pandas")
    fail_state = True

try:
    import pathlib
except Exception as e:
    print("Unable to import pathlib")
    fail_state = True

try:
    import re
except Exception as e:
    print("Unable to import re")
    fail_state = True

try:
    import sys
except Exception as e:
    print("Unable to import sys")
    fail_state = True

try:
    import tensorflow as tf
    print(tf.__version__)

except Exception as e:
    print("Unable to import tensorflow")
    fail_state = True

try:
    import random
except Exception as e:
    print("Unable to import random")
    fail_state = True

try:
    import datetime
except Exception as e:
    print("Unable to import datetime")
    fail_state = True

try:
    from matplotlib import pyplot as plt
except Exception as e:
    print("Unable to import matplotlib")
    fail_state = True

try:
    from PIL import Image
except Exception as e:
    print("Unable to import PIL")
    fail_state = True


if fail_state:
    print("Evn test FAILED")
else:
    print("Env Test complete: passed")
# if "tf" not in sys.modules:
#     print("TensorFlow not loaded")
# else:
#     print(tf.__version__)
