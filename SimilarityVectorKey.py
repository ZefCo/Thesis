import pandas
import numpy
import itertools
import pathlib


print_frame = pandas.DataFrame()


length = 9
rows = 1

weights = numpy.array([(1 / (2**(i))) for i in range(1, length + 1)])


insane = [numpy.reshape(numpy.array(i), (rows, length)) for i in itertools.product([0, 1], repeat = rows*length)]
insane_part2 = [numpy.sum(nuts * weights) for nuts in insane]
# print(len(insane_part2))
# print(insane_part2[2])
# print(insane[2][0])

insane = [thing[0] for thing in insane]
# print(test[2])
# print(type(insane))
# print(len(insane))
# print(numpy.sum(numpy.array(insane[0]) * weights))
# print(numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
# for a in range(0, 2):
#     for b in range(0, 2):
#         for c in range(0, 2):
#             for d in range(0, 2):
#                 for e in range(0, 2):
#                     for f in range(0, 2):
#                         for g in range(0, 2):
#                             for h in range(0, 2):
#                                 for i in range(0, 2):
#                                     for j in range(0, 2):
#                                         gdelta = numpy.array([a, b, c, d, e, f, g, h, i, j])
#                                         gstring = f'{a},{b},{c},{d},{e},{f},{g},{h},{i},{j}'
#                                         gnum = numpy.sum(weights * gdelta)
#                                         # all_values = numpy.concatenate((all_values, ))
#                                         # print(pandas.Series([gstring, values]))
#                                         # print(gnum, type(gnum))
#                                         # all_values = numpy.concatenate((all_values, values))
#                                         # print(pandas.Series([gstring, gnum]).T)
#                                         fuckthisgoddamnschool = pandas.Series([gnum, gstring])
#                                         print_frame = pandas.concat([print_frame, fuckthisgoddamnschool.to_frame().T], ignore_index = True)


# print(print_frame)
# rows, _ = print_frame.shape
# # print(rows)
# # print(print_frame)
# print_frame.to_csv(pathlib.Path.cwd() / 'Data_Files' / 'FUCKUH.csv', index=False)
# # print(count)
# # print(all_values)

print_frame = pandas.DataFrame(data = {"Similarity": insane, "G Score": insane_part2})
print_frame.to_csv(pathlib.Path.cwd() / f"Similarity Vector Key Length - {length}.csv", index = False)