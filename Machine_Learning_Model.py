import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.inspection import plot_partial_dependence
from sklearn.metrics import mean_squared_error
from keras.callbacks import EarlyStopping
from keras.callbacks import Callback
from keras import models
from keras import layers
from math import sqrt
import matplotlib.pyplot as plt
Import seaborn as sos


df = pd.read_csv('file_directory', low_memory = False)

y = df['Average '].astype('int64')

x = scale(x)
X_train, X_test, y_train, y_test = train_test_split(x, y, test_size = 0.2, random_state = 2)

model = models.Sequential()
parameters = {'activation':('relu', 'identity', 'logistic', 'tanh')}

model = GridSearchCV(model.add(layers.Dense(64, input_shape=(X_train.shape[1], )))
			model.add(layers.Dense(49))
			model.add(layers.Dense(1))
			model.compile(optimizer='Adadelta', loss='mean_squared_error',metrics=['accuracy'])
			, parameters)

model.summary()
history = model.fit(X_train, y_train, validation_split=0.25,epochs=150, batch_size=32, verbose=1)


from IPython.display import SVG
from keras.utils import model_to_dot

unseen = model.evaluate(X_test, y_test)
seen = model.evaluate(X_train, y_train)

print('Accuracy test data: ', round(unseen[1] * 100, 2),'%', '\n','Accuracy train data: ', round(seen[1] * 100, 2),'%')

hist = pd.DataFrame(history.history)

plt.plot(hist.axes[0].tolist(), hist['accuracy'])
plt.title('Model accuracy')
plt.ylabel('Accuracy')
plt.xlabel('Epoch')
plt.rcParams.update({'font.size': 25})

plt.show()
sns.pairplot(x[["logFC_TAM", "logFC_SUM", "logFC_MDALM2", "logFC_MDA"]], diag_kind="kde")
plt.show()
