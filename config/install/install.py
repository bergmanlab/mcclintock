URL = {
    "te-locate" : "https://downloads.sourceforge.net/project/te-locate/TE-locate.tar",
    "retroseq" : "https://github.com/tk2/RetroSeq/archive/9d4f3b5270af2383f40e6e7ea1204ea718365db2.zip",
    "temp" : "https://github.com/JialiUMassWengLab/TEMP/archive/4f67e1da836721a9f0999efa52e1e648fedb75fc.zip",
    "relocate" : "https://github.com/srobb1/RelocaTE/archive/ce3a2066e15f5c14e2887fdf8dce0485e1750e5b.zip",
    "ngs_te_mapper" : "https://github.com/bergmanlab/ngs_te_mapper/archive/fb23590200666fe66f1c417c5d5934385cb77ab9.zip",
    "popoolationte" : "http://downloads.sourceforge.net/project/popoolationte/popoolationte_1.02.zip",
    "popoolationte2": "http://downloads.sourceforge.net/project/popoolation-te2/popte2-v1.10.03.jar",
    "relocate2" : "https://github.com/stajichlab/RelocaTE2/archive/v2.0.1.zip"
}

MD5 = {
    "te-locate" : "c28e3dd7be89a8efcfbca8939923787c",
    "retroseq" : "54030657dd8476f1a65f5f71e809b4c3",
    "temp" : "f06e25244e9979510a9743bcd9c841a5",
    "relocate" : "64cc60f65154368518529fb57becc3d6",
    "ngs_te_mapper" : "971dccdb8112e754293f7b8b44799aeb",
    "popoolationte" : "315d266bd6da7487c919576dc46b92b4",
    "popoolationte2": "0c06186f0b1df215949f57a5cbb2d039",
    "relocate2" : "63eb7dbc09ffeb4063aef4857d198b06"
}


ENV_PATH = "{{envpath}}"
ENV = {
    "coverage" : ENV_PATH+"mcc_coverage.yml",
    "te-locate" : ENV_PATH+"mcc_telocate.yml",
    "retroseq" : ENV_PATH+"mcc_retroseq.yml",
    "temp" : ENV_PATH+"mcc_temp.yml",
    "relocate" : ENV_PATH+"mcc_relocate.yml",
    "ngs_te_mapper" : ENV_PATH+"mcc_ngs_te_mapper.yml",
    "popoolationte" : ENV_PATH+"mcc_popoolationte.yml",
    "popoolationte2" : ENV_PATH+"mcc_popoolationte2.yml",
    "relocate2" : ENV_PATH+"mcc_relocate2.yml",
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
    "popoolationte" : INSTALL_PATH+"tools/popoolationte/identify-te-insertsites.pl",
    "popoolationte2" : INSTALL_PATH+"tools/popoolationte2/popte2-v1.10.03.jar",
    "relocate2" : INSTALL_PATH+"tools/relocate2/relocaTE2.log",
    "processing" : INSTALL_PATH+"tools/processing.log"
}