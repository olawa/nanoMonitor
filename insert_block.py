                        try:
                            # Robust Chromosome Lookup
                            target_chrom = None
                            norm_chrom = chrom.replace("chr", "")
                            
                            if chrom in gene_models:
                                target_chrom = chrom
                            elif norm_chrom in gene_models:
                                target_chrom = norm_chrom
                            elif f"chr{norm_chrom}" in gene_models:
                                target_chrom = f"chr{norm_chrom}"
                            
                            if target_chrom:
                                overlaps = gene_models[target_chrom].overlap(start, end)
                                if overlaps:
                                    print(f"DEBUG: Overlaps for {chrom}:{start}-{end} -> {[o.data.get('name') for o in overlaps]}")
                                
                                best_gene = None
                                found_genes = set()
                                found_exons = defaultdict(set)
                                
                                for interval in overlaps:
                                    data = interval.data
                                    gname = data.get("name")
                                    ftype = data.get("type")
                                    if gname:
                                        found_genes.add(gname)
                                        best_gene = gname
                                        if ftype == "exon":
                                            enum = data.get("exon")
                                            if enum and enum != "?":
                                                found_exons[gname].add(enum)
                                
                                if found_genes:
                                    gene_strs = []
                                    for g in sorted(found_genes):
                                        if g in found_exons and found_exons[g]:
                                            try:
                                                exs = sorted([int(e) for e in found_exons[g]])
                                                if len(exs) > 1: ex_str = f"ex{min(exs)}-{max(exs)}"
                                                else: ex_str = f"ex{exs[0]}"
                                                gene_strs.append(f"{g}_{ex_str}")
                                            except:
                                                gene_strs.append(f"{g}_ex{','.join(sorted(found_exons[g]))}")
                                        else:
                                            gene_strs.append(g)
                                    gene_name = f"{chrom}:{start}-{end}({','.join(gene_strs)})"
                                elif best_gene:
                                    gene_name = best_gene
                        except Exception as e:
                            print(f"DEBUG: Gene resolution error: {e}")
                            pass
