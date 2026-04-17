import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns




control_file = "ENCFF622FSB.tsv"
kd_file = "ENCFF521KTC.tsv"

control = pd.read_csv(control_file, sep="\t")
kd = pd.read_csv(kd_file, sep="\t")
control = control[["gene_id", "expected_count"]]
kd = kd[["gene_id", "expected_count"]]
control.columns = ["gene_id", "control"]
kd.columns = ["gene_id", "knockdown"]


data = pd.merge(control, kd, on="gene_id")

data["control_log"] = np.log2(data["control"] + 1)
data["kd_log"] = np.log2(data["knockdown"] + 1)
data["log2FC"] = data["kd_log"] - data["control_log"]
data.to_csv("differential_expression_results.csv", index=False)

print("Results saved")




plt.figure(figsize=(8,6))
plt.scatter(data["log2FC"], data["kd_log"], alpha=0.3)

plt.xlabel("log2 Fold Change")
plt.ylabel("Expression (log2)")
plt.title("RNA-seq Expression Change")
plt.savefig("volcano_plot.png")
plt.close()



plt.figure(figsize=(6,6))
sns.boxplot(data=[data["control_log"], data["kd_log"]])

plt.xticks([0,1],["Control","TRA2A KD"])
plt.ylabel("log2 expression")

plt.title("Gene Expression Distribution")

plt.savefig("expression_boxplot.png")
plt.hist(data["control_log"],bins=50,alpha=0.5,label="Control")
plt.hist(data["kd_log"],bins=50,alpha=0.5,label="Knockdown")
plt.legend()
plt.title("QC Summary")
plt.savefig("qc_summary.png")

print("Plots saved")




top_genes = data.reindex(data["log2FC"].abs().sort_values(ascending=False).index)
top20 = top_genes.head(20)
top20.to_csv("top20_genes.csv", index=False)

print("Top 20 genes saved")




heatmap_data = top20[["control_log","kd_log"]]

plt.figure(figsize=(6,8))

sns.heatmap(
    heatmap_data,
    cmap="viridis",
    yticklabels=top20["gene_id"]
)

plt.title("Top 20 Differentially Expressed Genes")

plt.savefig("heatmap_top20_genes.png")

print("Heatmap saved")
