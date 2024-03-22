from matplotlib import pyplot as plt
import sys
import pandas
import pathlib
cwd = pathlib.Path.cwd()


log_file = sys.argv[1]
with open(log_file, "r") as lfile:
    for l in lfile:
        pass
    version_dir = pathlib.Path(l.rstrip())
    print(version_dir)

history_file = pathlib.Path(version_dir / "ModelSteps.csv")
history_data = pandas.read_csv(history_file, index_col="Unnamed: 0")

plt.plot(history_data["loss"])
plt.plot(history_data["val_loss"])
plt.xlabel("Epochs")
plt.legend(["Loss", "Valid Loss"])
plt.ylim(top = 10, bottom = 0)
plt.savefig(str(version_dir / "Loss_Graph_test.png"))

plt.close()

plt.plot(history_data["categorical_accuracy"])
plt.plot(history_data["val_categorical_accuracy"])
plt.xlabel("Epochs")
plt.legend(["Accuracy", "Valid Accuracy"])
plt.savefig(str(version_dir / "Accuracy_Graph_test.png"))