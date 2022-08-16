import numpy
# import pandas


def dissimilar_alignment(x, y, match = 0, mismatch = 1, gap = 2):
    # This was written with the idea of doing a G score in mind. Therefore this actually
    # looks for the MIN of the scores, penalizes gaps harshly, and adds a positive value to
    # a mismatch.
    
    # Establishing size of matrix
    row_x = len(x)
    row_y = len(y)

    # Creating initial matrix and adding in boundary conditions
    G = numpy.zeros((row_x + 1, row_y + 1))
    G[:, 0] = numpy.linspace(0, row_x * gap, row_x + 1)
    G[0, :] = numpy.linspace(0, row_y * gap, row_y + 1)
    # print(G)

    # Creating pointers to trace through alignment
    P = numpy.zeros((row_x + 1, row_y + 1))
    # setting first column to point up, first row to point to the left
    P[:, 0] = 3
    P[0, :] = 4

    # Traceback Scores
    t = numpy.zeros(3)
    for i in range(row_x):
        for j in range(row_y):
            # establish if it's a match or a mismatch
            if x[i] == y[j]:
                t[0] = G[i, j] + match
            else:
                t[0] = G[i, j] + mismatch
            
            # Establish the gap scores
            t[1] = G[i, j + 1] + gap
            t[2] = G[i + 1, j] + gap
            
            # Take the higest of the scores
            tmin = numpy.min(t)
            
            # Put the min in that position in G
            G[i + 1, j + 1] = tmin
            # print(f"G matrix:\n{G}")

            if t[0] == tmin:
                # Arrow points diagonally
                P[i + 1, j + 1] += 2
            if t[1] == tmin:
                # Arrow points to up
                P[i + 1, j + 1] += 3
            if t[2] == tmin:
                # Arrow points to the left
                P[i + 1, j + 1] += 4
            # print(f"{P}")

    # print(G[row_x, row_y])
    # Trace through optimal alignment
    i, j = row_x, row_y
    rx, ry = [], []

    while i > 0 or j > 0:
        # 2, 3, 4, 5, 6, 7, 9 represent all possible values of 2, 3, and 4 put together
        # This checks to see which one it is, with preference to the diagonal arrows
        if P[i, j] in [2, 5, 6, 9]:
            # print(P[i, j])
            rx.append(x[i - 1])
            ry.append(y[j - 1])
            i -= 1
            j -= 1

        elif P[i, j] in [3, 5, 7, 9]:
            # print(P[i, j])
            rx.append(x[i - 1])
            ry.append('-')
            i -= 1
        elif P[i, j] in [4, 6, 7, 9]:
            # print(P[i, j])
            rx.append('-')
            ry.append(y[j - 1])
            j -= 1
        
    rx = ''.join(rx)[::-1]
    ry = ''.join(ry)[::-1]

    return rx, ry


if __name__ in '__main__':
    # x = 'TTATTGTACC'
    # y = 'TAGTTTGTGA'
    y = "CGTTACGTAA"
    x = "CGTAC"
    # y = 'TCCCAG'

    sequence = dissimilar_alignment(x, y)
    # print(sequence)
    print(f"{sequence[0]}\n{sequence[1]}")
