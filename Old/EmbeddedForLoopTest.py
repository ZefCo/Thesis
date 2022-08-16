from numpy import isin


tries = 10
denominators = 10

for d in range(denominators):
    
    for t in range(tries):
        try:
            a = 10 / d
            break
        except ZeroDivisionError as e:
            error_message = e
            print(f"Failed to find a, trying again, attempt: {t + 1}")
        except Exception as e:
            error_message = e
            print(f"Writing Error {error_message} type {type(error_message)} to log")
            break
    else:
        print(f"Writing {error_message} to log")
        continue


    if isinstance(a, float):
        if a % 2 != 0:
            continue
        else:
            print(f"{d} clears all conditions: a = {a}")    
