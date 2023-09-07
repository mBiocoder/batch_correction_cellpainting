library(data.table)
library(ggplot2)
library(funkyheatmap)
library(dplyr)
library(tibble)

setwd("../../data/scib_metrics_results/")

palettes <- tribble(
  ~palette, ~colours,
  "Blues", grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues") %>% c("#011636")))(101),
  "Reds", grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Reds")[-8:-9]))(101)
)
column_groups <- tribble(
  ~Category, ~group, ~palette,
  "Batch Correction", "group1", "Blues",
  "Biological Conservation", "group2", "Reds"
) 



# direct integration
dt <- fread("jump_scib_direct.csv")
dt <- melt(dt, id.vars = c("Metric"), variable.name="Method", value.name="Value")
dt <- dcast(dt, Method~Metric, value.var="Value")
colnames(dt)[1] <- "id"

dt <- rbind(dt, list("min", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
dt <- rbind(dt, list("max", 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))

row_info <- tribble(
  ~id,
  "Unintegrated",
  "Harmony",
  "Scanorama",
  "scGen",
  "scVI",
  "scANVI",
  "gaushVI",
  "gaushANVI"
)

column_info <- tribble(
  ~id,  ~group, ~name, ~geom, ~palette, ~legend,
  "id", "", "", "text", NA, FALSE,
  "PCR_comparison", "group1", "PCR batch", "funkyrect", "Blues", TRUE,
  "iLISI", "group1", "Graph iLISI", "funkyrect", "Blues", FALSE,
  "Graph_connectivity", "group1",  "Graph connectivity", "funkyrect", "Blues", FALSE,
  "Batch_ASW", "group1", "Batch ASW", "funkyrect", "Blues", FALSE,
  "kBET", "group1", "kBET", "funkyrect", "Blues", FALSE,
  "NMI", "group2", "NMI cluster/label", "funkyrect", "Reds", FALSE,
  "ARI", "group2", "ARI cluster/label", "funkyrect", "Reds", FALSE,
  "Isolated_labels_F1", "group2", "Isolated label F1", "funkyrect", "Reds", FALSE,
  "Silhouette", "group2", "Isolated label silhouette", "funkyrect", "Reds", FALSE,
  "Isolated_labels_ASW", "group2", "Cell type ASW", "funkyrect", "Reds", FALSE,
  "cLISI", "group2", "Graph cLISI", "funkyrect", "Reds", FALSE
)
funky_heatmap(dt, column_info = column_info, palettes = palettes,
              column_groups = column_groups, row_info = row_info)

ggsave("scib_direct_test.png", bg="white", dpi=400, units = "in",
       width=9, height=9)




# high-level integration
dt <- fread("jump_scib_high.csv")
dt <- melt(dt, id.vars = c("Metric"), variable.name="Method", value.name="Value")
dt <- dcast(dt, Method~Metric, value.var="Value")
colnames(dt)[1] <- "id"

dt <- rbind(dt, list("min", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
dt <- rbind(dt, list("max", 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))

row_info <- tribble(
  ~id,
  "Unintegrated",
  "Harmony",
  "Scanorama",
  "scGen",
  "scVI",
  "scANVI",
  "gaushANVI"
)

column_info <- tribble(
  ~id,  ~group, ~name, ~geom, ~palette, ~legend,
  "id", "", "", "text", NA, FALSE,
  "PCR_comparison", "group1", "PCR batch", "funkyrect", "Blues", TRUE,
  "iLISI", "group1", "Graph iLISI", "funkyrect", "Blues", FALSE,
  "Graph_connectivity", "group1",  "Graph connectivity", "funkyrect", "Blues", FALSE,
  "Batch_ASW", "group1", "Batch ASW", "funkyrect", "Blues", FALSE,
  "kBET", "group1", "kBET", "funkyrect", "Blues", FALSE,
  "NMI", "group2", "NMI cluster/label", "funkyrect", "Reds", FALSE,
  "ARI", "group2", "ARI cluster/label", "funkyrect", "Reds", FALSE,
  "Isolated_labels_F1", "group2", "Isolated label F1", "funkyrect", "Reds", FALSE,
  "Silhouette", "group2", "Isolated label silhouette", "funkyrect", "Reds", FALSE,
  "Isolated_labels_ASW", "group2", "Cell type ASW", "funkyrect", "Reds", FALSE,
  "cLISI", "group2", "Graph cLISI", "funkyrect", "Reds", FALSE
)
funky_heatmap(dt, column_info = column_info, palettes = palettes,
                   column_groups = column_groups, row_info = row_info)

ggsave("scib_high.png", bg="white", dpi="print", units = "in",
       width=8, height=7.3)




# low-level integration
dt <- fread("jump_scib_low.csv")
dt <- melt(dt, id.vars = c("Metric", "Source"), variable.name="Method", value.name="Value")
dt <- dcast(dt, Method+Source~Metric, value.var="Value")
dt[, id:=rownames(dt)]
dt[, Source:=factor(dt$Source, levels=c("source_2", "source_3", "source_4",
                                        "source_5", "source_6", "source_7",
                                        "source_8", "source_9", "source_10",
                                        "source_11", "source_13"))]

dt <- rbind(dt, list("min", "", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ""))
dt <- rbind(dt, list("max", "", 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ""))

d <- dt[, c("id", "Source", "Method")]
d <- d[Source!="",]
d[, Method:=factor(d$Method, levels=c("Unintegrated", "Harmony", "Scanorama", "scGen", "scVI", "scANVI"))]
setorder(d, cols = "Source", "Method")
d[, group:="group2"]
d[Source=="source_3", group:="group3"]
d[Source=="source_4", group:="group4"]
d[Source=="source_5", group:="group5"]
d[Source=="source_6", group:="group6"]
d[Source=="source_7", group:="group7"]
d[Source=="source_8", group:="group8"]
d[Source=="source_9", group:="group9"]
d[Source=="source_10", group:="group10"]
d[Source=="source_11", group:="group11"]
d[Source=="source_13", group:="group13"]
row_info <- d[, c("id", "group")]

row_groups <- tribble(
  ~group, ~level1,
  "group2", "Source 2",
  "group3", "Source 3",
  "group4", "Source 4",
  "group5", "Source 5",
  "group6", "Source 6",
  "group7", "Source 7",
  "group8", "Source 8",
  "group9", "Source 9",
  "group10", "Source 10",
  "group11", "Source 11",
  "group13", "Source 13"
) 

column_info <- tribble(
  ~id,  ~group, ~name, ~geom, ~palette, ~legend,
  "Method", "", "", "text", NA, FALSE,
  "PCR_comparison", "group1", "PCR batch", "funkyrect", "Blues", TRUE,
  "iLISI", "group1", "Graph iLISI", "funkyrect", "Blues", FALSE,
  "Graph_connectivity", "group1",  "Graph connectivity", "funkyrect", "Blues", FALSE,
  "Batch_ASW", "group1", "Batch ASW", "funkyrect", "Blues", FALSE,
  "kBET", "group1", "kBET", "funkyrect", "Blues", FALSE,
  "NMI", "group2", "NMI cluster/label", "funkyrect", "Reds", FALSE,
  "ARI", "group2", "ARI cluster/label", "funkyrect", "Reds", FALSE,
  "Isolated_labels_F1", "group2", "Isolated label F1", "funkyrect", "Reds", FALSE,
  "Silhouette", "group2", "Isolated label silhouette", "funkyrect", "Reds", FALSE,
  "Isolated_labels_ASW", "group2", "Cell type ASW", "funkyrect", "Reds", FALSE,
  "cLISI", "group2", "Graph cLISI", "funkyrect", "Reds", FALSE
)
funky_heatmap(dt, column_info = column_info, palettes = palettes,
                   column_groups = column_groups, row_info = row_info,
                   row_groups = row_groups)

ggsave("scib_low.png", bg="white", dpi="print", units = "in",
       width=15, height=25)


funky_heatmap(dt, column_info = column_info, palettes = palettes,
              column_groups = column_groups, row_info = row_info[group=="group13",],
              row_groups = data.table(row_groups)[group=="group13",])
ggsave("scib_low_13.png", bg="white", dpi="print", units = "in",
       width=9.17, height=6.77)



# low-level integration, mean over all sources
mean_dt <- dt[ , list(ARI=mean(ARI), Batch_ASW=mean(Batch_ASW),
           Graph_connectivity=mean(Graph_connectivity),
           Isolated_labels_ASW=mean(Isolated_labels_ASW),
           Isolated_labels_F1=mean(Isolated_labels_F1),
           NMI=mean(NMI), PCR_comparison=mean(PCR_comparison),
           Silhouette=mean(Silhouette), cLISI=mean(cLISI),
           iLISI=mean(iLISI), kBET=mean(kBET)), by="Method"]
colnames(mean_dt)[1] <- "id"
row_info <- tribble(
  ~id,
  "Unintegrated",
  "Harmony",
  "Scanorama",
  "scGen",
  "scVI",
  "scANVI"
)

column_info <- tribble(
  ~id,  ~group, ~name, ~geom, ~palette, ~legend,
  "id", "", "", "text", NA, FALSE,
  "PCR_comparison", "group1", "PCR batch", "funkyrect", "Blues", TRUE,
  "iLISI", "group1", "Graph iLISI", "funkyrect", "Blues", FALSE,
  "Graph_connectivity", "group1",  "Graph connectivity", "funkyrect", "Blues", FALSE,
  "Batch_ASW", "group1", "Batch ASW", "funkyrect", "Blues", FALSE,
  "kBET", "group1", "kBET", "funkyrect", "Blues", FALSE,
  "NMI", "group2", "NMI cluster/label", "funkyrect", "Reds", FALSE,
  "ARI", "group2", "ARI cluster/label", "funkyrect", "Reds", FALSE,
  "Isolated_labels_F1", "group2", "Isolated label F1", "funkyrect", "Reds", FALSE,
  "Silhouette", "group2", "Isolated label silhouette", "funkyrect", "Reds", FALSE,
  "Isolated_labels_ASW", "group2", "Cell type ASW", "funkyrect", "Reds", FALSE,
  "cLISI", "group2", "Graph cLISI", "funkyrect", "Reds", FALSE
)
funky_heatmap(mean_dt, column_info = column_info, palettes = palettes,
                   column_groups = column_groups, row_info = row_info)

ggsave("scib_low_mean.png", bg="white", dpi="print", units = "in",
       width=9.17, height=6.77)
ggsave("scib_low_mean_test.png", bg="white", dpi=400, units = "in",
       width=9, height=9)







# direct integration, low-level metrics
dt <- fread("jump_scib_direct_low.csv")
dt <- melt(dt, id.vars = c("Metric", "Source"), variable.name="Method", value.name="Value")
dt <- dcast(dt, Method+Source~Metric, value.var="Value")
dt[, id:=rownames(dt)]
dt[, Source:=factor(dt$Source, levels=c("source_2", "source_3", "source_4",
                                        "source_5", "source_6", "source_7",
                                        "source_8", "source_9", "source_10",
                                        "source_11", "source_13"))]

dt <- rbind(dt, list("min", "", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ""))
dt <- rbind(dt, list("max", "", 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ""))

d <- dt[, c("id", "Source", "Method")]
d <- d[Source!="",]
d[, Method:=factor(d$Method, levels=c("Unintegrated", "Harmony", "Scanorama", "scGen", "scVI", "scANVI"))]
setorder(d, cols = "Source", "Method")
d[, group:="group2"]
d[Source=="source_3", group:="group3"]
d[Source=="source_4", group:="group4"]
d[Source=="source_5", group:="group5"]
d[Source=="source_6", group:="group6"]
d[Source=="source_7", group:="group7"]
d[Source=="source_8", group:="group8"]
d[Source=="source_9", group:="group9"]
d[Source=="source_10", group:="group10"]
d[Source=="source_11", group:="group11"]
d[Source=="source_13", group:="group13"]
row_info <- d[, c("id", "group")]

row_groups <- tribble(
  ~group, ~level1,
  "group2", "Source 2",
  "group3", "Source 3",
  "group4", "Source 4",
  "group5", "Source 5",
  "group6", "Source 6",
  "group7", "Source 7",
  "group8", "Source 8",
  "group9", "Source 9",
  "group10", "Source 10",
  "group11", "Source 11",
  "group13", "Source 13"
) 

column_info <- tribble(
  ~id,  ~group, ~name, ~geom, ~palette, ~legend,
  "Method", "", "", "text", NA, FALSE,
  "PCR_comparison", "group1", "PCR batch", "funkyrect", "Blues", TRUE,
  "iLISI", "group1", "Graph iLISI", "funkyrect", "Blues", FALSE,
  "Graph_connectivity", "group1",  "Graph connectivity", "funkyrect", "Blues", FALSE,
  "Batch_ASW", "group1", "Batch ASW", "funkyrect", "Blues", FALSE,
  "kBET", "group1", "kBET", "funkyrect", "Blues", FALSE,
  "NMI", "group2", "NMI cluster/label", "funkyrect", "Reds", FALSE,
  "ARI", "group2", "ARI cluster/label", "funkyrect", "Reds", FALSE,
  "Isolated_labels_F1", "group2", "Isolated label F1", "funkyrect", "Reds", FALSE,
  "Silhouette", "group2", "Isolated label silhouette", "funkyrect", "Reds", FALSE,
  "Isolated_labels_ASW", "group2", "Cell type ASW", "funkyrect", "Reds", FALSE,
  "cLISI", "group2", "Graph cLISI", "funkyrect", "Reds", FALSE
)
g <- funky_heatmap(dt, column_info = column_info, palettes = palettes,
                   column_groups = column_groups, row_info = row_info,
                   row_groups = row_groups)
g
#ggsave("scib_direct_low.png", g, bg="white", dpi="print", units = "in",
#       width=15, height=25)


funky_heatmap(dt, column_info = column_info, palettes = palettes,
              column_groups = column_groups, row_info = row_info[group=="group13",],
              row_groups = data.table(row_groups)[group=="group13",])
ggsave("scib_direct_low_13.png", g, bg="white", dpi="print", units = "in",
       width=9.17, height=6.77)



# low-level integration, mean over all sources
mean_dt <- dt[ , list(ARI=mean(ARI), Batch_ASW=mean(Batch_ASW),
                      Graph_connectivity=mean(Graph_connectivity),
                      Isolated_labels_ASW=mean(Isolated_labels_ASW),
                      Isolated_labels_F1=mean(Isolated_labels_F1),
                      NMI=mean(NMI), PCR_comparison=mean(PCR_comparison),
                      Silhouette=mean(Silhouette), cLISI=mean(cLISI),
                      iLISI=mean(iLISI), kBET=mean(kBET)), by="Method"]
colnames(mean_dt)[1] <- "id"
row_info <- tribble(
  ~id,
  "Unintegrated",
  "Harmony",
  "Scanorama",
  "scGen",
  "scVI",
  "scANVI"
)

column_info <- tribble(
  ~id,  ~group, ~name, ~geom, ~palette, ~legend,
  "id", "", "", "text", NA, FALSE,
  "PCR_comparison", "group1", "PCR batch", "funkyrect", "Blues", TRUE,
  "iLISI", "group1", "Graph iLISI", "funkyrect", "Blues", FALSE,
  "Graph_connectivity", "group1",  "Graph connectivity", "funkyrect", "Blues", FALSE,
  "Batch_ASW", "group1", "Batch ASW", "funkyrect", "Blues", FALSE,
  "kBET", "group1", "kBET", "funkyrect", "Blues", FALSE,
  "NMI", "group2", "NMI cluster/label", "funkyrect", "Reds", FALSE,
  "ARI", "group2", "ARI cluster/label", "funkyrect", "Reds", FALSE,
  "Isolated_labels_F1", "group2", "Isolated label F1", "funkyrect", "Reds", FALSE,
  "Silhouette", "group2", "Isolated label silhouette", "funkyrect", "Reds", FALSE,
  "Isolated_labels_ASW", "group2", "Cell type ASW", "funkyrect", "Reds", FALSE,
  "cLISI", "group2", "Graph cLISI", "funkyrect", "Reds", FALSE
)
funky_heatmap(mean_dt, column_info = column_info, palettes = palettes,
                   column_groups = column_groups, row_info = row_info)

ggsave("scib_direct_low_mean.png", g, bg="white", dpi="print", units = "in",
       width=9.17, height=6.77)







### Own method
### JUMP
# correct for source
dt <- fread("scib_own_jump_source.csv")
dt <- melt(dt, id.vars = c("Metric"), variable.name="Method", value.name="Value")
dt <- dcast(dt, Method~Metric, value.var="Value")
colnames(dt)[1] <- "id"

dt <- rbind(dt, list("min", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
dt <- rbind(dt, list("max", 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))

row_info <- tribble(
  ~id,
  "Unintegrated",
  "gaushVI",
  "gaushANVI"
)

column_info <- tribble(
  ~id,  ~group, ~name, ~geom, ~palette, ~legend,
  "id", "", "", "text", NA, FALSE,
  "PCR_comparison", "group1", "PCR batch", "funkyrect", "Blues", TRUE,
  "iLISI", "group1", "Graph iLISI", "funkyrect", "Blues", FALSE,
  "Graph_connectivity", "group1",  "Graph connectivity", "funkyrect", "Blues", FALSE,
  "Batch_ASW", "group1", "Batch ASW", "funkyrect", "Blues", FALSE,
  "kBET", "group1", "kBET", "funkyrect", "Blues", FALSE,
  "NMI", "group2", "NMI cluster/label", "funkyrect", "Reds", FALSE,
  "ARI", "group2", "ARI cluster/label", "funkyrect", "Reds", FALSE,
  "Isolated_labels_F1", "group2", "Isolated label F1", "funkyrect", "Reds", FALSE,
  "Silhouette", "group2", "Isolated label silhouette", "funkyrect", "Reds", FALSE,
  "Isolated_labels_ASW", "group2", "Cell type ASW", "funkyrect", "Reds", FALSE,
  "cLISI", "group2", "Graph cLISI", "funkyrect", "Reds", FALSE
)
funky_heatmap(dt, column_info = column_info, palettes = palettes,
                   column_groups = column_groups, row_info = row_info)
ggsave("scib_jump_source.png", bg="white", dpi="print", units = "in",
       width=9.17, height=6.77)





# RxRx19b
dt <- fread("scib_own_rxrx19b.csv")
dt <- melt(dt, id.vars = c("Metric"), variable.name="Method", value.name="Value")
dt <- dcast(dt, Method~Metric, value.var="Value")
colnames(dt)[1] <- "id"

dt <- rbind(dt, list("min", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
dt <- rbind(dt, list("max", 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))

row_info <- tribble(
  ~id,
  "Unintegrated",
  "Harmony",
  "Scanorama",
  "scGen",
  "scVI",
  "scANVI",
  "gaushVI",
  "gaushANVI"
)

column_info <- tribble(
  ~id,  ~group, ~name, ~geom, ~palette, ~legend,
  "id", "", "", "text", NA, FALSE,
  "iLISI", "group1", "Graph iLISI", "funkyrect", "Blues", FALSE,
  "Graph_connectivity", "group1",  "Graph connectivity", "funkyrect", "Blues", FALSE,
  "Batch_ASW", "group1", "Batch ASW", "funkyrect", "Blues", FALSE,
  "NMI", "group2", "NMI cluster/label", "funkyrect", "Reds", FALSE,
  "ARI", "group2", "ARI cluster/label", "funkyrect", "Reds", FALSE,
  "Isolated_labels_F1", "group2", "Isolated label F1", "funkyrect", "Reds", FALSE,
  "Silhouette", "group2", "Isolated label silhouette", "funkyrect", "Reds", FALSE,
  "Isolated_labels_ASW", "group2", "Cell type ASW", "funkyrect", "Reds", FALSE,
  "cLISI", "group2", "Graph cLISI", "funkyrect", "Reds", FALSE
)
column_groups <- tribble(
  ~Category, ~group, ~palette,
  "Batch Corr.", "group1", "Blues",
  "Bio Conservation", "group2", "Reds"
) 
funky_heatmap(dt, column_info = column_info, palettes = palettes,
              column_groups = column_groups, row_info = row_info)

ggsave("own_metrics/scib_rxrx19b.png", bg="white", dpi=400, units = "in",
       width=8, height=7.5)




# RxRx1
dt <- fread("own_metrics/scib_rxrx1_controls.csv")
dt <- melt(dt, id.vars = c("Metric"), variable.name="Method", value.name="Value")
dt <- dcast(dt, Method~Metric, value.var="Value")
colnames(dt)[1] <- "id"

dt <- rbind(dt, list("min", 0, 0, 0, 0, 0, 0, 0, 0, 0))
dt <- rbind(dt, list("max", 1, 1, 1, 1, 1, 1, 1, 1, 1))

row_info <- tribble(
  ~id,
  "Unintegrated",
  "Harmony",
  "Scanorama",
  "scGen",
  "scVI",
  "scANVI",
  "gaushVI",
  "gaushANVI"
)

column_info <- tribble(
  ~id,  ~group, ~name, ~geom, ~palette, ~legend,
  "id", "", "", "text", NA, FALSE,
  "PCR_comparison", "group1", "PCR batch", "funkyrect", "Blues", TRUE,
  "iLISI", "group1", "Graph iLISI", "funkyrect", "Blues", FALSE,
  "Graph_connectivity", "group1",  "Graph connectivity", "funkyrect", "Blues", FALSE,
  "Batch_ASW", "group1", "Batch ASW", "funkyrect", "Blues", FALSE,
  #"kBET", "group1", "kBET", "funkyrect", "Blues", FALSE,
  "NMI", "group2", "NMI cluster/label", "funkyrect", "Reds", FALSE,
  "Silhouette", "group2", "Isolated label silhouette", "funkyrect", "Reds", FALSE,
  "ARI", "group2", "ARI cluster/label", "funkyrect", "Reds", FALSE,
  "cLISI", "group2", "Graph cLISI", "funkyrect", "Reds", FALSE
)
column_groups <- tribble(
  ~Category, ~group, ~palette,
  "Batch Correction", "group1", "Blues",
  "Bio Conservation", "group2", "Reds"
) 
funky_heatmap(dt, column_info = column_info, palettes = palettes,
              column_groups = column_groups, row_info = row_info)

ggsave("own_metrics/scib_rxrx1_controls.png", bg="white", dpi=400, units = "in",
       width=8, height=7.3)




### BBBC021
dt <- fread("own_metrics/scib_bbbc021.csv")
dt <- melt(dt, id.vars = c("Metric"), variable.name="Method", value.name="Value")
dt <- dcast(dt, Method~Metric, value.var="Value")
colnames(dt)[1] <- "id"

dt <- rbind(dt, list("min", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
dt <- rbind(dt, list("max", 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))

row_info <- tribble(
  ~id,
  "Unintegrated",
  "scVI",
  "scANVI"
)

column_info <- tribble(
  ~id,  ~group, ~name, ~geom, ~palette, ~legend,
  "id", "", "", "text", NA, FALSE,
  "PCR_comparison", "group1", "PCR batch", "funkyrect", "Blues", TRUE,
  "iLISI", "group1", "Graph iLISI", "funkyrect", "Blues", FALSE,
  "Graph_connectivity", "group1",  "Graph connectivity", "funkyrect", "Blues", FALSE,
  "Batch_ASW", "group1", "Batch ASW", "funkyrect", "Blues", FALSE,
  "kBET", "group1", "kBET", "funkyrect", "Blues", FALSE,
  "NMI", "group2", "NMI cluster/label", "funkyrect", "Reds", FALSE,
  "ARI", "group2", "ARI cluster/label", "funkyrect", "Reds", FALSE,
  "Isolated_labels_F1", "group2", "Isolated label F1", "funkyrect", "Reds", FALSE,
  "Silhouette", "group2", "Isolated label silhouette", "funkyrect", "Reds", FALSE,
  "Isolated_labels_ASW", "group2", "Cell type ASW", "funkyrect", "Reds", FALSE,
  "cLISI", "group2", "Graph cLISI", "funkyrect", "Reds", FALSE
)
column_groups <- tribble(
  ~Category, ~group, ~palette,
  "Batch Correction", "group1", "Blues",
  "Biological Conservation", "group2", "Reds"
) 
funky_heatmap(dt, column_info = column_info, palettes = palettes,
              column_groups = column_groups, row_info = row_info)

ggsave("own_metrics/scib_bbbc021.png", bg="white", dpi=400, units = "in",
       width=8, height=7.3)



### JUMP directly integrate for plate
dt <- fread("own_metrics/scib_own_jump_plate.csv")
dt <- melt(dt, id.vars = c("Metric"), variable.name="Method", value.name="Value")
dt <- dcast(dt, Method~Metric, value.var="Value")
colnames(dt)[1] <- "id"

dt <- rbind(dt, list("min", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
dt <- rbind(dt, list("max", 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))

row_info <- tribble(
  ~id,
  "Unintegrated",
  "gaushVI",
  "gaushANVI"
)

column_info <- tribble(
  ~id,  ~group, ~name, ~geom, ~palette, ~legend,
  "id", "", "", "text", NA, FALSE,
  "PCR_comparison", "group1", "PCR batch", "funkyrect", "Blues", TRUE,
  "iLISI", "group1", "Graph iLISI", "funkyrect", "Blues", FALSE,
  "Graph_connectivity", "group1",  "Graph connectivity", "funkyrect", "Blues", FALSE,
  "Batch_ASW", "group1", "Batch ASW", "funkyrect", "Blues", FALSE,
  "kBET", "group1", "kBET", "funkyrect", "Blues", FALSE,
  "NMI", "group2", "NMI cluster/label", "funkyrect", "Reds", FALSE,
  "ARI", "group2", "ARI cluster/label", "funkyrect", "Reds", FALSE,
  "Isolated_labels_F1", "group2", "Isolated label F1", "funkyrect", "Reds", FALSE,
  "Silhouette", "group2", "Isolated label silhouette", "funkyrect", "Reds", FALSE,
  "Isolated_labels_ASW", "group2", "Cell type ASW", "funkyrect", "Reds", FALSE,
  "cLISI", "group2", "Graph cLISI", "funkyrect", "Reds", FALSE
)
column_groups <- tribble(
  ~Category, ~group, ~palette,
  "Batch Correction", "group1", "Blues",
  "Biological Conservation", "group2", "Reds"
) 
funky_heatmap(dt, column_info = column_info, palettes = palettes,
              column_groups = column_groups, row_info = row_info)

ggsave("own_metrics/scib_jump_plate.png", bg="white", dpi=400, units = "in",
       width=9, height=9)
