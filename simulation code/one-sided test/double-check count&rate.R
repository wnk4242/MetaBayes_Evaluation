

#这个Script检查False/True Success之类的count数目是否正确
#rowise数据集里包含了FS/TS等metrics已经数好的数目
#比如FEMABF这个method，在以下simulation条件下的false success的个数为284
rowisecount <- MABF_rates_null_rowise_BFcutoff1 %>% 
  filter(method == "FEMABF", orig.alpha == 0.05, orig.n == 20, bias.level == "low", rep.number == 2, rep.n == 40)
sum(rowisecount$FS2)

#我们再从raw data里数一遍
#这是与上面的数据集对应的raw data
rawBF <- FEMABF_lists_0.2null_regrouped_deltap[["0,0.2_20_none_low"]][["2_40"]]
#我们选取所有前500 row的数据，这500 row是根据underlying effect size is 0生成的
#column 2是original study的p值，column 3是original study的effect size，剩下的是replication studies的Baye factor值
#因为要数多少false success的数目，我们只选符合p <.05 and es > 0这种条件的row
raw_pesBF <- rawBF[1:500, 2:503] %>% 
  as_tibble() %>% 
  filter(original_p<0.05 & observed_es > 0)
#选好这些row后，数有多少BF>1，这些BF符合FS2的标准
raw.BF <- raw_pesBF[,4:502]
#sum(rowisecount$FS2)和sum(raw.BF>1)得到的结果都是284个，说明我们的数据是正确的
sum(raw.BF>1)

#也可以用下面这个思路得到同样的结果
count_FS2 <- FEMABF_lists_0.2null_regrouped_deltap[["0,0.2_20_none_low"]][["2_40"]][1:500, ] %>%
  as_tibble() %>%
  filter(original_p < 0.05, observed_es > 0) %>%
  select(4:ncol(.)) %>%
  summarise(total = sum(across(everything(), ~ .x > 1))) %>%
  pull(total)
