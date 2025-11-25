library(dplyr)
library(stringr)


pattern <- Sys.getenv('pattern') 
path_out <- Sys.getenv('path')

#LOAD ALL OUTPUT
setwd(path_out)
files <- list.files(pattern=paste0(pattern, ".Rda"))
print(files)

out_all <- data.frame()
pars_all <- data.frame()
ystart_all <- data.frame()
yend_all <- data.frame()
for(i in 1:length(files)) {
  file <- files[i]
  print(file)
  load(file)
  out_all <- bind_rows(out_all, as.data.frame(post$out) %>% mutate(chain=i))
  pars_all <- bind_rows(pars_all, as.data.frame(post$resample) %>% mutate(chain=i))
  ystart_all <- bind_rows(ystart_all, as.data.frame(post$ystart_all) %>% mutate(chain=i))
  yend_all <- bind_rows(yend_all, as.data.frame(post$yend_all) %>% mutate(chain=i))
}

write.csv(pars_all, paste0("params_post_", pattern, ".csv"), row.names=F)
write.csv(out_all, paste0("out_post_", pattern, ".csv"), row.names=F)
write.csv(ystart_all, paste0("ystart_post_", pattern, ".csv"), row.names=F)
write.csv(yend_all, paste0("yend_post_", pattern, ".csv"), row.names=F)