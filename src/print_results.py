
from mynotebook import get_joined, show_all_orfs

# import matplotlib as mpl
# mpl.use("Agg")  # Use a backend that does not support on-screen

joined = get_joined("cfeintact/plasma")

#
# Size cutoffs are determined manually.
#

# print("###########")
# print("## Sizes ##")
# print("###########")
# show_all_orfs(joined, "CFEIntact", "intact", "size", 0.01)

print("###############")
print("## Distances ##")
print("###############")
show_all_orfs(joined, "CFEIntact", "intact", "distance", 0.0001)

print("##################")
print("## Indel impact ##")
print("##################")
show_all_orfs(joined, "CFEIntact", "intact", "indel impact", 0.0001)
