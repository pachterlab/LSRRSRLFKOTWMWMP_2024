library(bambu)
mOC1 <- "/home/rebekah/LRGASP_data/mouse_ES_ONT_CapTrap_1_IQ/mouse_ES_ONT_CapTrap_1_IQ/aux/mouse_ES_ONT_CapTrap_1_IQ_ENCFF356OJC_9a93d3_4c0465_9384c8.bam"
mOC2 <- "/home/rebekah/LRGASP_data/mouse_ES_ONT_CapTrap_2_IQ/mouse_ES_ONT_CapTrap_2_IQ/aux/mouse_ES_ONT_CapTrap_2_IQ_ENCFF275RMO_2456ac_0d2b46_3988ae.bam"
mOC3 <- "/home/rebekah/LRGASP_data/mouse_ES_ONT_CapTrap_3_IQ/mouse_ES_ONT_CapTrap_3_IQ/aux/mouse_ES_ONT_CapTrap_3_IQ_ENCFF056EOI_3f7b8c_ccc3fe_0d006e.bam"
mOc1 <- "/home/rebekah/LRGASP_data/mouse_ES_ONT_cDNA_1_IQ/mouse_ES_ONT_cDNA_1_IQ/aux/mouse_ES_ONT_cDNA_1_IQ_ENCFF683TBO_b510f9_62acea_f4b355.bam"
mOc2 <- "/home/rebekah/LRGASP_data/mouse_ES_ONT_cDNA_2_IQ/mouse_ES_ONT_cDNA_2_IQ/aux/mouse_ES_ONT_cDNA_2_IQ_ENCFF232YSU_99e16c_ad191d_e94380.bam"
mOc3 <- "/home/rebekah/LRGASP_data/mouse_ES_ONT_cDNA_3_IQ/mouse_ES_ONT_cDNA_3_IQ/aux/mouse_ES_ONT_cDNA_3_IQ_ENCFF288PBL_cde748_d99675_94011b.bam"
mOd1 <- "/home/rebekah/LRGASP_data/mouse_ES_ONT_dRNA_1_IQ/mouse_ES_ONT_dRNA_1_IQ/aux/mouse_ES_ONT_dRNA_1_IQ_ENCFF765AEC.perlm_2f5c20_a5fb53_f9b00b.bam"
mOd2 <- "/home/rebekah/LRGASP_data/mouse_ES_ONT_dRNA_2_IQ/mouse_ES_ONT_dRNA_2_IQ/aux/mouse_ES_ONT_dRNA_2_IQ_ENCFF349BIN.perlm_e52efa_2b279c_90496c.bam"
mOd3 <- "/home/rebekah/LRGASP_data/mouse_ES_ONT_dRNA_3_IQ/mouse_ES_ONT_dRNA_3_IQ/aux/mouse_ES_ONT_dRNA_3_IQ_ENCFF412NKJ.perlm_24fcdb_18da91_12cec6.bam"
mPC1 <- "/home/rebekah/LRGASP_data/mouse_ES_PacBio_CapTrap_1_IQ/mouse_ES_PacBio_CapTrap_1_IQ/aux/mouse_ES_PacBio_CapTrap_1_IQ_ENCFF535DQA_268127_17b0cb_c3282e.bam"
mPC2 <- "/home/rebekah/LRGASP_data/mouse_ES_PacBio_CapTrap_2_IQ/mouse_ES_PacBio_CapTrap_2_IQ/aux/mouse_ES_PacBio_CapTrap_2_IQ_ENCFF310IPO_9706ff_5f5f58_3e26d0.bam"
mPC3 <- "/home/rebekah/LRGASP_data/mouse_ES_PacBio_CapTrap_3_IQ/mouse_ES_PacBio_CapTrap_3_IQ/aux/mouse_ES_PacBio_CapTrap_3_IQ_ENCFF654JHQ_9bbafd_4a8ea3_800793.bam"
mPc1 <- "/home/rebekah/LRGASP_data/mouse_ES_PacBio_cDNA_1_IQ/mouse_ES_PacBio_cDNA_1_IQ/aux/mouse_ES_PacBio_cDNA_1_IQ_ENCFF874VSI_211a0b_4f6dce_1a48da.bam"
mPc2 <- "/home/rebekah/LRGASP_data/mouse_ES_PacBio_cDNA_2_IQ/mouse_ES_PacBio_cDNA_2_IQ/aux/mouse_ES_PacBio_cDNA_2_IQ_ENCFF667VXS_643fba_4ef602_16bd11.bam"
mPc3 <- "/home/rebekah/LRGASP_data/mouse_ES_PacBio_cDNA_3_IQ/mouse_ES_PacBio_cDNA_3_IQ/aux/mouse_ES_PacBio_cDNA_3_IQ_ENCFF313VYZ_1fcc8d_853b20_976e76.bam"
fa.file <- "/home/rebekah/lrgasp_grcm39_sirvs.fasta"

gtf.file <- "/home/rebekah/lrgasp_gencode_vM27_sirvs.gtf"

bambuAnnotations <- prepareAnnotations(gtf.file)

se <- bambu(reads = c(mOc1, mOc2, mOc3), annotations = bambuAnnotations, genome = fa.file, discovery = FALSE)
writeBambuOutput(se, path = "LRGASP_data/er_quant_tcc/mouse_ES_ONT_cDNA_bambu")

se <- bambu(reads = c(mOC1, mOC2, mOC3), annotations = bambuAnnotations, genome = fa.file, discovery = FALSE)
writeBambuOutput(se, path = "LRGASP_data/er_quant_tcc/mouse_ES_ONT_CapTrap_bambu")

se <- bambu(reads = c(mOd1, mOd2, mOd3), annotations = bambuAnnotations, genome = fa.file, discovery = FALSE)
writeBambuOutput(se, path = "LRGASP_data/er_quant_tcc/mouse_ES_ONT_dRNA_bambu")

se <- bambu(reads = c(mPC1, mPC2, mPC3), annotations = bambuAnnotations, genome = fa.file, discovery = FALSE)
writeBambuOutput(se, path = "LRGASP_data/er_quant_tcc/mouse_ES_PacBio_CapTrap_bambu")

se <- bambu(reads = c(mPc1, mPc2, mPc3), annotations = bambuAnnotations, genome = fa.file, discovery = FALSE)
writeBambuOutput(se, path = "LRGASP_data/er_quant_tcc/mouse_ES_PacBio_cDNA_bambu")

hOC1 <- "/home/rebekah/LRGASP_data/human_WTC11_ONT_CapTrap_1_IQ/human_WTC11_ONT_CapTrap_1_IQ/aux/human_WTC11_ONT_CapTrap_1_IQ_ENCFF654SNK_7804a4_075ecf_0a496a.bam"
hOC2 <- "/home/rebekah/LRGASP_data/human_WTC11_ONT_CapTrap_2_IQ/human_WTC11_ONT_CapTrap_2_IQ/aux/human_WTC11_ONT_CapTrap_2_IQ_ENCFF934KDM_238b5c_75b0ab_58987c.bam"
hOC3 <- "/home/rebekah/LRGASP_data/human_WTC11_ONT_CapTrap_3_IQ/human_WTC11_ONT_CapTrap_3_IQ/aux/human_WTC11_ONT_CapTrap_3_IQ_ENCFF104BNW_a19dad_276fd5_055076.bam"
hOc1 <- "/home/rebekah/LRGASP_data/human_WTC11_ONT_cDNA_1_IQ/human_WTC11_ONT_cDNA_1_IQ/aux/human_WTC11_ONT_cDNA_1_IQ_ENCFF263YFG_eee77d_e805e6_900730.bam"
hOc2 <- "/home/rebekah/LRGASP_data/human_WTC11_ONT_cDNA_2_IQ/human_WTC11_ONT_cDNA_2_IQ/aux/human_WTC11_ONT_cDNA_2_IQ_ENCFF023EXJ_1efed2_3ee9b5_8050b5.bam"
hOc3 <- "/home/rebekah/LRGASP_data/human_WTC11_ONT_cDNA_3_IQ/human_WTC11_ONT_cDNA_3_IQ/aux/human_WTC11_ONT_cDNA_3_IQ_ENCFF961HLO_b487b8_bbf723_b50b2f.bam"
hOd1 <- "/home/rebekah/LRGASP_data/human_WTC11_ONT_dRNA_1_IQ/human_WTC11_ONT_dRNA_1_IQ/aux/human_WTC11_ONT_dRNA_1_IQ_ENCFF155CFF.perlm_a16a78_ed7b19_4f0e56.bam"
hOd2 <- "/home/rebekah/LRGASP_data/human_WTC11_ONT_dRNA_2_IQ/human_WTC11_ONT_dRNA_2_IQ/aux/human_WTC11_ONT_dRNA_2_IQ_ENCFF771DIX.perlm_a82dce_be598f_c1f000.bam"
hOd3 <- "/home/rebekah/LRGASP_data/human_WTC11_ONT_dRNA_3_IQ/human_WTC11_ONT_dRNA_3_IQ/aux/human_WTC11_ONT_dRNA_3_IQ_ENCFF600LIU.perlm_6c73a6_09971f_4a0f51.bam"
hPC1 <- "/home/rebekah/LRGASP_data/human_WTC11_PacBio_CapTrap_1_IQ/human_WTC11_PacBio_CapTrap_1_IQ/aux/human_WTC11_PacBio_CapTrap_1_IQ_ENCFF105WIJ_55f08a_316a57_16b28d.bam"
hPC2 <- "/home/rebekah/LRGASP_data/human_WTC11_PacBio_CapTrap_2_IQ/human_WTC11_PacBio_CapTrap_2_IQ/aux/human_WTC11_PacBio_CapTrap_2_IQ_ENCFF212HLP_118ec9_9e93b2_031bb1.bam"
hPC3 <- "/home/rebekah/LRGASP_data/human_WTC11_PacBio_CapTrap_3_IQ/human_WTC11_PacBio_CapTrap_3_IQ/aux/human_WTC11_PacBio_CapTrap_3_IQ_ENCFF003QZT_43fd19_6899dd_1c25a0.bam"
hPc1 <- "/home/rebekah/LRGASP_data/human_WTC11_PacBio_cDNA_1_IQ/human_WTC11_PacBio_cDNA_1_IQ/aux/human_WTC11_PacBio_cDNA_1_IQ_ENCFF563QZR_96aec9_e708b5_e37015.bam"
hPc2 <- "/home/rebekah/LRGASP_data/human_WTC11_PacBio_cDNA_2_IQ/human_WTC11_PacBio_cDNA_2_IQ/aux/human_WTC11_PacBio_cDNA_2_IQ_ENCFF370NFS_e02a76_006501_fa67bd.bam"
hPc3 <- "/home/rebekah/LRGASP_data/human_WTC11_PacBio_cDNA_3_IQ/human_WTC11_PacBio_cDNA_3_IQ/aux/human_WTC11_PacBio_cDNA_3_IQ_ENCFF245IPA_2dbff5_49a37a_fd68b6.bam"

fa.file <- "/home/rebekah/lrgasp_grch38_sirvs.fasta"

gtf.file <- "/home/rebekah/lrgasp_gencode_v38_sirvs.gtf"

bambuAnnotations <- prepareAnnotations(gtf.file)

se <- bambu(reads = c(hOd1, hOd2, hOd3), annotations = bambuAnnotations, genome = fa.file, discovery = FALSE)
writeBambuOutput(se, path = "LRGASP_data/er_quant_tcc/human_WTC11_ONT_dRNA_bambu")

se <- bambu(reads = c(hPC1, hPC2, hPC3), annotations = bambuAnnotations, genome = fa.file, discovery = FALSE)
writeBambuOutput(se, path = "LRGASP_data/er_quant_tcc/human_WTC11_PacBio_CapTrap_bambu")

se <- bambu(reads = c(hPc1, hPc2, hPc3), annotations = bambuAnnotations, genome = fa.file, discovery = FALSE)
writeBambuOutput(se, path = "LRGASP_data/er_quant_tcc/human_WTC11_PacBio_cDNA_bambu")

se <- bambu(reads = c(hOC1, hOC2, hOC3), annotations = bambuAnnotations, genome = fa.file, discovery = FALSE)
writeBambuOutput(se, path = "LRGASP_data/er_quant_tcc/human_WTC11_ONT_CapTrap_bambu")

se <- bambu(reads = c(hOc1, hOc2, hOc3), annotations = bambuAnnotations, genome = fa.file, discovery = FALSE)
writeBambuOutput(se, path = "LRGASP_data/er_quant_tcc/human_WTC11_ONT_cDNA_bambu")
