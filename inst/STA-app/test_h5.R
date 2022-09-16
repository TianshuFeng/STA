library(rhdf5)

obj_mapper <- readRDS("/Volumes/GoogleDrive/My Drive/Mayo_intern files/Projects/Mapper_Visualization/Package/SemiMapper/inst/STA-app/chick_foot.rds")
descript <- read.csv("chick_foot.csv")
tp_data = chicken_generator(1)

test_network <- readRDS("test_network.RDS")
tp_data <- read.table("normcount/normcount.txt")
save_network_h5(obj_mapper = test_network,
                dataset = as.data.frame(t(tp_data)),
                file = "test.h5")


save_network_h5(obj_mapper = obj_mapper,
                      dataset = tp_data,
                      file = "test.h5")

a <- load_network_h5(file = "test.h5")

aaa <- read.table("conditions.txt", stringsAsFactors = F)

table(aaa[test_network$points_in_vertex[[44]],2])

h5createFile("test.h5")

h5ls("test.h5")


obj_mapper_list <- obj_mapper
attr(obj_mapper_list, "class") <- NULL

h5write(obj_mapper_list, file = "test.h5", name = "test")
a <- h5read(file = "test.h5", name = "test")
a$points_in_level_set <- a$points_in_level_set[mixedsort(names(a$points_in_level_set))]

ii <- sapply(descript, is.factor)
descript[ii] <- lapply(descript[ii], as.character)

h5closeAll()
h5delete(file = "test.h5", name = "descript")
h5write(descript, file = "test.h5", name = "descript",DataFrameAsCompound = FALSE)

h5ls("test.h5")
b <- h5read(file = "test.h5", name = "descript/Group")


h5write("", file = "test.h5", name = "testnull")

aa <- data.frame(name = names(descript),
                 numeric = sapply(descript, is.numeric),
                 stringsAsFactors = FALSE)
h5closeAll()
h5delete(file = "test.h5", name = "colname_type1")
h5write(aa, file = "test.h5", name = "colname_type1")
bb <- h5read(file = "test.h5", name = "dataset")
bb

h5f = H5Fopen("test.h5")
h5f&'descript'&'Group'[]



dom_grp <- c()
for (i in obj_mapper$points_in_vertex) {
  dom_grp <-
    c(dom_grp, names(sort(table(b[i]), decreasing = T))[1])
}
