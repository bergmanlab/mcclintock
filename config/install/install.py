URL = {
    "te-locate" : "https://downloads.sourceforge.net/project/te-locate/TE-locate.tar",
    "retroseq" : "https://github.com/tk2/RetroSeq/archive/9d4f3b5270af2383f40e6e7ea1204ea718365db2.zip",
    "temp" : "https://github.com/JialiUMassWengLab/TEMP/archive/4f67e1da836721a9f0999efa52e1e648fedb75fc.zip",
    "relocate" : "https://github.com/srobb1/RelocaTE/archive/ce3a2066e15f5c14e2887fdf8dce0485e1750e5b.zip",
    "ngs_te_mapper" : "https://github.com/bergmanlab/ngs_te_mapper/archive/f9f48996ac346ac86d57edbd00534aa1227b753e.zip",
    "ngs_te_mapper2": "https://github.com/bergmanlab/ngs_te_mapper/archive/98fad0c92e2d3089d96dc85967bfbcef917922af.zip",
    "popoolationte" : "http://downloads.sourceforge.net/project/popoolationte/popoolationte_1.02.zip",
    "popoolationte2": "http://downloads.sourceforge.net/project/popoolation-te2/popte2-v1.10.03.jar",
    "relocate2" : "https://github.com/stajichlab/RelocaTE2/archive/v2.0.1.zip",
    "tepid" : "https://github.com/ListerLab/TEPID/archive/ad46d65b5c41bf8a9171215d49b3ffaecdceaab0.zip",
    "teflon": "https://github.com/jradrion/TEFLoN/archive/9eca0152f3dd9dc6c44787a30d590f7e321b442c.zip",
    "jitterbug": "https://github.com/elzbth/jitterbug/archive/b6b3f9c7ee4af042d4410137269f174a1399b752.zip"
}

MD5 = {
    "te-locate" : "c28e3dd7be89a8efcfbca8939923787c",
    "retroseq" : "54030657dd8476f1a65f5f71e809b4c3",
    "temp" : "f06e25244e9979510a9743bcd9c841a5",
    "relocate" : "64cc60f65154368518529fb57becc3d6",
    "ngs_te_mapper" : "a00e9962ee3daca9fac915c23c6af203",
    "ngs_te_mapper2": "18db9da34f1728db86e65f67e9fe646c",
    "popoolationte" : "315d266bd6da7487c919576dc46b92b4",
    "popoolationte2": "0c06186f0b1df215949f57a5cbb2d039",
    "relocate2" : "63eb7dbc09ffeb4063aef4857d198b06",
    "tepid" : "51f61acf138b07842074f054aba780d5",
    "teflon" : "87080f4f455a97eaeaa5d8e3d8c4f35f",
    "jitterbug": "05f165fb6a79a7a611d09981a355d27e"
}


ENV_PATH = "{{envpath}}"
ENV = {
    "coverage" : ENV_PATH+"mcc_coverage.yml",
    "te-locate" : ENV_PATH+"mcc_telocate.yml",
    "retroseq" : ENV_PATH+"mcc_retroseq.yml",
    "temp" : ENV_PATH+"mcc_temp.yml",
    "relocate" : ENV_PATH+"mcc_relocate.yml",
    "ngs_te_mapper" : ENV_PATH+"mcc_ngs_te_mapper.yml",
    "ngs_te_mapper2": ENV_PATH+"mcc_ngs_te_mapper2.yml",
    "popoolationte" : ENV_PATH+"mcc_popoolationte.yml",
    "popoolationte2" : ENV_PATH+"mcc_popoolationte2.yml",
    "relocate2" : ENV_PATH+"mcc_relocate2.yml",
    "tepid" : ENV_PATH+"mcc_tepid.yml",
    "teflon" : ENV_PATH+"mcc_teflon.yml",
    "jitterbug": ENV_PATH+"mcc_jitterbug.yml",
    "setup_reads": ENV_PATH+"mcc_setup_reads.yml",
    "mcc_processing":  ENV_PATH+"mcc_processing.yml"
}


INSTALL_PATH = "{{inspath}}"
OUTPUT = {
    "coverage" : INSTALL_PATH+"tools/coverage/coverage.log",
    "te-locate" : INSTALL_PATH+"tools/te-locate/TE_locate.pl",
    "retroseq" : INSTALL_PATH+"tools/retroseq/bin/retroseq.pl",
    "temp" : INSTALL_PATH+"tools/temp/scripts/TEMP_Insertion.sh",
    "relocate" : INSTALL_PATH+"tools/relocate/scripts/relocaTE_insertionFinder.pl",
    "ngs_te_mapper" : INSTALL_PATH+"tools/ngs_te_mapper/sourceCode/run_ngs_te_mapper.sh",
    "ngs_te_mapper2": INSTALL_PATH+"tools/ngs_te_mapper2/sourceCode/ngs_te_mapper.py",
    "popoolationte" : INSTALL_PATH+"tools/popoolationte/identify-te-insertsites.pl",
    "popoolationte2" : INSTALL_PATH+"tools/popoolationte2/popte2-v1.10.03.jar",
    "relocate2" : INSTALL_PATH+"tools/relocate2/relocaTE2.log",
    "tepid" : INSTALL_PATH+"tools/tepid/tepid-map",
    "teflon" : INSTALL_PATH+"tools/teflon/teflon.v0.4.py",
    "jitterbug": INSTALL_PATH+"tools/jitterbug/jitterbug.py",
    "processing" : INSTALL_PATH+"tools/processing.log",
    "setup_reads" : INSTALL_PATH+"tools/setup_reads.log"
}