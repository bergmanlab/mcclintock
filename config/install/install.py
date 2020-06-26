URL = {
    "te-locate" : "https://downloads.sourceforge.net/project/te-locate/TE-locate.tar",
    "retroseq" : "https://github.com/tk2/RetroSeq/archive/700d4f76a3b996686652866f2b81fefc6f0241e0.zip",
    "temp" : "https://github.com/JialiUMassWengLab/TEMP/archive/d2500b904e2020d6a1075347b398525ede5feae1.zip",
    "relocate" : "https://github.com/srobb1/RelocaTE/archive/ce3a2066e15f5c14e2887fdf8dce0485e1750e5b.zip",
    "ngs_te_mapper" : "https://github.com/bergmanlab/ngs_te_mapper/archive/fb23590200666fe66f1c417c5d5934385cb77ab9.zip",
    "popoolationte" : "http://downloads.sourceforge.net/project/popoolationte/popoolationte_1.02.zip",
    "popoolationte2": "http://downloads.sourceforge.net/project/popoolation-te2/popte2-v1.10.03.jar",
    "relocate2" : "https://github.com/stajichlab/RelocaTE2/archive/v2.0.1.zip"
}

MD5 = {
    "te-locate" : "c28e3dd7be89a8efcfbca8939923787c",
    "retroseq" : "6141e5aa4285125ab1145d1db53c4114",
    "temp" : "cfb77b2d1e82e38971e793f5e71f7d2c",
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