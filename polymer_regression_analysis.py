import lime
import grakel
import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn import linear_model
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
import sklearn.ensemble
import pyGPs


# Truncates strings up to the character char.
def cleanData(df, char, exclude_columns):
    for col in df:
        if col in exclude_columns:
            #  We may not want to truncate values in some columns, such as block ratios
            continue
        if col == 'Monomers' or col == 'Big_Smile' or col == 'Block_ratio' or col == 'Paper':
            continue
        for number in range(len(df[col])):
            finalValue = ''
            if isinstance(df[col][number], float) or isinstance(df[col][number], int):
                df[col][number] = float(df[col][number])
            else:
                for c in df[col][number]:
                    if c in char:
                        break
                    else:
                        if c != '>': # to replace entries such as '>1800' with '1800'
                            finalValue += c
                df[col][number] = float(finalValue)
    return df

def cleanDataWithVariances(df, char, exclude_columns):
    df_variances = df.copy()
    for col in df:
        if col in exclude_columns:
            #  We may not want to truncate values in some columns, such as block ratios
            continue
        if col == 'Monomers' or col == 'Big_Smile' or col == 'Block_ratio' or col == 'Paper':
            continue
        for number in range(len(df[col])):
            finalValue = ''
            finalValueVariance = ''
            # If entry is a float already, we record it and assign its variance to zero:
            if isinstance(df[col][number], float) or isinstance(df[col][number], int):
                df[col][number] = float(df[col][number])
                if not np.isnan(df[col][number]):
                    df_variances[col][number] = 0
            else:
                is_var = False
                for c in df[col][number]:
                    if c in char:
                        if c != '±': # We truncate the number as it is.
                            break
                        else: # We star recording the variance from next character
                            is_var = True
                          
                    else:
                        if not is_var:
                            if c != '>': # to replace entries such as '>1800' with '1800'
                                finalValue += c
                        else: 
                            finalValueVariance += c
                if not is_var: # If value doesn't have a variance, assign zero to it.
                    finalValueVariance = '0'
                    
                            
                df[col][number] = float(finalValue)
#                 if is_var:
                df_variances[col][number] = float(finalValueVariance)
    return df, df_variances

def normalise(df):
    for col in df:
        if col == 'Monomers' or col == 'Big_Smile' or col == 'Block_ratio' or col == 'Paper':
            continue
        # Create x, where x the 'scores' column's values as floats
        x = df[[col]].values.astype(float)

        # Create a minimum and maximum processor object
        min_max_scaler = preprocessing.MinMaxScaler()

        # Create an object to transform the data to fit minmax processor
        x_scaled = min_max_scaler.fit_transform(x)

        # Run the normalizer on the dataframe
        df[col] = pd.DataFrame(x_scaled)


def normalise_standard(df):
    for col in df:
        if col == 'Monomers' or col == 'Big Smile' or col == 'Block ratio' or col == 'Paper':
            continue
        # Create x, where x the 'scores' column's values as floats
        x = df[[col]].values.astype(float)

        # Create a minimum and maximum processor object
        stdard_scaler = preprocessing.StandardScaler()

        # Create an object to transform the data to fit standard processor
        x_scaled = stdard_scaler.fit_transform(x)

        # Run the normalizer on the dataframe
        df[col] = pd.DataFrame(x_scaled)


def calculateErrorInMethod(X, y, method = 'rf', error= 'mean_sq', runs = 1000):
    mean_sq_error =[]
    for run in range(runs):
        X_train, X_test, y_train, y_test = train_test_split(X, y,  test_size=0.3, random_state=run)
        if method == 'rf':
            met = sklearn.ensemble.RandomForestRegressor(n_estimators=1000)
            met.fit(X_train, y_train)
            y_pred = met.predict(X_test)
        elif method == 'gp':
            #  Make sure we are  not using sparse matrices, or even dense matrices. We need np arrays for GP to work.
            if type(X_train) != np.ndarray:
                x = X_train.toarray()
            else:
                x = X_train
            if type(X_test) != np.ndarray:
                z = X_test.toarray()
            else:
                z = X_test
            model = pyGPs.GPR()
            m = pyGPs.mean.Linear( D=x.shape[1] )
            k = pyGPs.cov.RBF()
            y_train = np.array(y_train)
            model.getPosterior(x, y_train)
#             model.setPrior(mean=m, kernel=k)
#             print('mean is', model.meanfunc.getMean())
            model.setPrior(mean=m)
            model.setOptimizer("Minimize", num_restarts=50)
            model.optimize()
            y_pred = model.predict(z)[0]
        else:
            met = linear_model.LinearRegression()
            met.fit(X_train, y_train)
            y_pred = met.predict(X_test)
        if error == 'mean_sq':
            if runs == 1:
                print("Real:", y_test)
                print("Pred:", np.around(y_pred, decimals = 1))
            mean_sq_error.append(mean_squared_error(y_test, y_pred))
        if error == 'r2':
            mean_sq_error.append(r2_score(y_test, y_pred))
            print(y_test, y_pred)
    #     print(np.mean(mean_sq_error))
    #     print(np.std(mean_sq_error))
    return np.mean(mean_sq_error), np.std(mean_sq_error)




