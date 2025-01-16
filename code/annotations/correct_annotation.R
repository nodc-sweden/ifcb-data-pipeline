library(iRfcb)

ifcb_py_install(envname = ".virtualenvs/iRfcb")

ifcb_path <- Sys.getenv("ifcb_path")

ifcb_correct_annotation(file.path(ifcb_path, "manual", "Tångesund"),
                        file.path(ifcb_path, "manual", "Tångesund"),
                        file.path("output", "Scrippsiella_group_selected_images_scrippsiella_round.txt"),
                        44)
