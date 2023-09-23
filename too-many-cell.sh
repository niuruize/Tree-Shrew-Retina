################################################################################
# too-many-cell analysis
################################################################################
%%R
# too-many-cell-counts
snRNA <- subset(snRNA, downsample=1000)
count = as.data.frame(as.matrix(snRNA@assays$integrated@data))
write.csv(count,"count.csv")
label_celltype_species <- data.frame(item=colnames(count), label=snRNA$celltype_species)
write.csv(label_celltype_species,"labels.csv", row.names = F)

%%bash
docker run -it --rm -v /Users/niuruize/Downloads/TOO:/TOO \
      gregoryschwartz/too-many-cells:2.0.0.0 make-tree \
    --matrix-path /TOO/input/count.csv \
    --labels-file /TOO/input/labels.csv \
    --smart-cutoff 1 \
    --min-size 1 \
    --draw-collection "PieRing" \
    --dendrogram-output "dendrogram.pdf" \
    --output /TOO/out \
    > clusters.csv
