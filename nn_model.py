import joblib
import numpy as np
import pandas as pd
from keras.backend import sigmoid
from keras.layers import Dense
from keras.models import Sequential
from tensorflow.keras.optimizers import Adam
from keras.utils.generic_utils import get_custom_objects
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

df = pd.read_csv("training1.csv", header=0)
df2 = pd.read_csv("training2.csv", header=0)
df3 = pd.read_csv("training3.csv", header=0)
df4 = pd.read_csv("training4.csv", header=0)
df5 = pd.read_csv("training5.csv", header=0)
df6 = pd.read_csv("training6.csv", header=0)

X = df4.iloc[:, 0:19].values
y = df4.iloc[:, 19].values
X2 = df5.iloc[:, 0:19].values
y2 = df5.iloc[:, 19].values

X = np.concatenate((X, X2), axis=0)
y = np.concatenate((y, y2), axis=0)

train_data = X
train_targets = y
train_data = train_data.astype(float)

scaler = StandardScaler()
scaler.fit(train_data)
scaled_X = scaler.transform(train_data)

partial_train_data, val_data, partial_tran_targets, val_targets = train_test_split(scaled_X, train_targets,
                                                                                   test_size=0.2)

X1 = partial_train_data
Y1 = partial_tran_targets

X2 = val_data
Y2 = val_targets

print('len', len(X1), len(X2))


def swish(x):
    return x * sigmoid(1 * x)


def create_model(X1, Y1, X2, Y2):
    neurons = 256
    get_custom_objects().update({'swish': swish})
    model = Sequential()
    model.add(Dense(512, input_shape=[19], activation="linear"))
    model.add(Dense(neurons, activation="swish"))
    model.add(Dense(neurons, activation="swish"))
    model.add(Dense(32, activation="swish"))
    model.add(Dense(1, activation="linear"))
    model.compile(loss='mean_squared_error', optimizer=Adam(learning_rate=1e-3, decay=1e-3 / 200), metrics=["accuracy"])
    history = model.fit(X1, Y1, validation_data=(X2, Y2), epochs=500, batch_size=64, verbose=True)
    return model, history, neurons


scaler_filename = "scaler_jobswish800ref.save"

joblib.dump(scaler, scaler_filename)
model, history, neurons = create_model(X1, Y1, X2, Y2)
model.save("my_modelswish800ref")
PredValSet = model.predict(X2)
# print(PredValSet)
# print(np.linalg.norm(Y2, PredValSet))

plt.xlabel('epoch')
plt.ylabel('loss_magnitude')
plt.plot(history.history['loss'])
plt.grid('minor')
plt.savefig(rf'myfig_{neurons}_neurons.png')
# plt.show()
