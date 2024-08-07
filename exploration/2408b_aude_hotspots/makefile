joint_graph_fld := results/paeruginosa_ST235/backbone_joints/joints_pangraph
tree := results/paeruginosa_ST235/pangraph/coretree.nwk
jpos := results/paeruginosa_ST235/backbone_joints/joints_pos.json
gbks := data/gbk
ann := results/paeruginosa_ST235/annotations/junct_pos
gen_lengths := results/paeruginosa_ST235/pangraph/genome_lengths.csv


data:
	mkdir -p data

data/tree.nwk: data
	ln -sf ../../../$(tree) data/tree.nwk

data/joints_pos.json: data
	ln -sf ../../../$(jpos) data/joints_pos.json

data/gbk: data
	ln -sf ../../../$(gbks) data/gbk

data/genome_lengths.csv: data
	ln -sf ../../../$(gen_lengths) data/genome_lengths.csv

data/ann: data
	mkdir -p data/ann

data/ann/%.csv: data/ann
	ln -sf ../../../../$(ann)/$*_real.csv data/ann/$*.csv

res/%/joint_graph.json:
	mkdir -p res/$*
	ln -sf ../../../../$(joint_graph_fld)/$*.json $@

res/%/gbk_annotations.csv res/%/tool_annotations.csv: data/genome_lengths.csv data/gbk data/ann/defensefinder.csv data/ann/genomad.csv data/ann/integronfinder.csv data/ann/ISEScan.csv data/joints_pos.json res/%/joint_graph.json
	python3 01_extract_annotations.py --hs $*

res/%/annotations.pdf: res/%/gbk_annotations.csv res/%/tool_annotations.csv res/%/joint_graph.json data/tree.nwk
	python3 02_draw_all_genes.py --hs $*

LAB := \
	EBPKATYDVF_r__HFYDQJUCFM_f \
	MVIUGRLZNX_r__QHNTPCOKJC_r \
	

TGTS := $(foreach lab,$(LAB),res/$(lab)/annotations.pdf)

.SECONDARY:

.PHONY: all
all:
	$(MAKE) $(TGTS)