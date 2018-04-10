mat = snakemake@input[["mat_ip"]]
output = snakemake@output[[1]]
control = snakemake@params[["control"]]

print(mat)
print(control)
print(output)
print(class(control))

#write.table(print(class(control)),output)