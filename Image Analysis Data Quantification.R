#Install and load packages

install.packages("tidyverse")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("patchwork")

library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)

#Import Cell Profiler data

Control_Nuclei=read.csv("C:/Users/grant/OneDrive - Nexus365/Cowley Lab/Columbus Data/Control Nuclei.csv")
Day_2_Nuclei=read.csv("C:/Users/grant/OneDrive - Nexus365/Cowley Lab/Columbus Data/Day 2 Nuclei.csv")
Day_4_Nuclei=read.csv("C:/Users/grant/OneDrive - Nexus365/Cowley Lab/Columbus Data/Day 4 Nuclei.csv")
Day_7_Nuclei=read.csv("C:/Users/grant/OneDrive - Nexus365/Cowley Lab/Columbus Data/Day 7 Nuclei.csv")
Day_14_Nuclei=read.csv("C:/Users/grant/OneDrive - Nexus365/Cowley Lab/Columbus Data/Day 14 Nuclei.csv")

#Clean data to include the Well ID and Mean PU1 Signal Intensity of each object (nucleus)

Control_Nuclei_1=Control_Nuclei%>%
  select(Metadata_Well,Intensity_MeanIntensity_PU1)%>%
  mutate(Plate=1)%>%
  rename(Well=Metadata_Well,PU1_MeanIntensity=Intensity_MeanIntensity_PU1)

Day_2_Nuclei_1=Day_2_Nuclei%>%
  select(Metadata_Well,Intensity_MeanIntensity_PU1)%>%
  mutate(Plate=2)%>%
  rename(Well=Metadata_Well,PU1_MeanIntensity=Intensity_MeanIntensity_PU1)

Day_4_Nuclei_1=Day_4_Nuclei%>%
  select(Metadata_Well,Intensity_MeanIntensity_PU1)%>%
  mutate(Plate=3)%>%
  rename(Well=Metadata_Well,PU1_MeanIntensity=Intensity_MeanIntensity_PU1)

Day_7_Nuclei_1=Day_7_Nuclei%>%
  select(Metadata_Well,Intensity_MeanIntensity_PU1)%>%
  mutate(Plate=4)%>%
  rename(Well=Metadata_Well,PU1_MeanIntensity=Intensity_MeanIntensity_PU1)

Day_14_Nuclei_1=Day_14_Nuclei%>%
  select(Metadata_Well,Intensity_MeanIntensity_PU1)%>%
  mutate(Plate=5)%>%
  rename(Well=Metadata_Well,PU1_MeanIntensity=Intensity_MeanIntensity_PU1)

#Calculate signal intensity cutoffs for PU1 positivity using 2ndry Ab-only control wells
#Cutoffs will be 2*SD of mean PU1 signal intensity in nuclei stained with 2ndry Ab only
#There will be a separate 2ndry Ab-only control well for each cell line on each plate

Control_Nuclei_2ndry_Only_Stats=Control_Nuclei_1%>%
  filter(Well=="B003"|Well=="C003")%>%
  summarise(PU1_Mean = mean(PU1_MeanIntensity),PU1_SD = sd(PU1_MeanIntensity))%>%
  mutate(PU1_Cutoff = PU1_Mean + 2 * PU1_SD)
Control_Nuclei_PU1_Cutoff=Control_Nuclei_2ndry_Only_Stats$PU1_Cutoff  

Day_2_Nuclei_2ndry_Only_Stats_LP=Day_2_Nuclei_1%>%
  filter(Well=="C003")%>%
  summarise(PU1_Mean = mean(PU1_MeanIntensity),PU1_SD = sd(PU1_MeanIntensity))%>%
  mutate(PU1_Cutoff = PU1_Mean + 2 * PU1_SD)
Day_2_Nuclei_LP_PU1_Cutoff=Day_2_Nuclei_2ndry_Only_Stats_LP$PU1_Cutoff

Day_2_Nuclei_2ndry_Only_Stats_KOLF2.1S=Day_2_Nuclei_1%>%
  filter(Well=="C007")%>%
  summarise(PU1_Mean = mean(PU1_MeanIntensity),PU1_SD = sd(PU1_MeanIntensity))%>%
  mutate(PU1_Cutoff = PU1_Mean + 2 * PU1_SD)
Day_2_Nuclei_KOLF2.1S_PU1_Cutoff=Day_2_Nuclei_2ndry_Only_Stats_KOLF2.1S$PU1_Cutoff

Day_4_Nuclei_2ndry_Only_Stats_LP=Day_4_Nuclei_1%>%
  filter(Well=="C003")%>%
  summarise(PU1_Mean = mean(PU1_MeanIntensity),PU1_SD = sd(PU1_MeanIntensity))%>%
  mutate(PU1_Cutoff = PU1_Mean + 2 * PU1_SD)
Day_4_Nuclei_LP_PU1_Cutoff=Day_4_Nuclei_2ndry_Only_Stats_LP$PU1_Cutoff

Day_4_Nuclei_2ndry_Only_Stats_KOLF2.1S=Day_4_Nuclei_1%>%
  filter(Well=="C007")%>%
  summarise(PU1_Mean = mean(PU1_MeanIntensity),PU1_SD = sd(PU1_MeanIntensity))%>%
  mutate(PU1_Cutoff = PU1_Mean + 2 * PU1_SD)
Day_4_Nuclei_KOLF2.1S_PU1_Cutoff=Day_4_Nuclei_2ndry_Only_Stats_KOLF2.1S$PU1_Cutoff

Day_7_Nuclei_2ndry_Only_Stats_LP=Day_7_Nuclei_1%>%
  filter(Well=="C003")%>%
  summarise(PU1_Mean = mean(PU1_MeanIntensity),PU1_SD = sd(PU1_MeanIntensity))%>%
  mutate(PU1_Cutoff = PU1_Mean + 2 * PU1_SD)
Day_7_Nuclei_LP_PU1_Cutoff=Day_7_Nuclei_2ndry_Only_Stats_LP$PU1_Cutoff

Day_7_Nuclei_2ndry_Only_Stats_KOLF2.1S=Day_7_Nuclei_1%>%
  filter(Well=="C007")%>%
  summarise(PU1_Mean = mean(PU1_MeanIntensity),PU1_SD = sd(PU1_MeanIntensity))%>%
  mutate(PU1_Cutoff = PU1_Mean + 2 * PU1_SD)
Day_7_Nuclei_KOLF2.1S_PU1_Cutoff=Day_7_Nuclei_2ndry_Only_Stats_KOLF2.1S$PU1_Cutoff

Day_14_Nuclei_2ndry_Only_Stats_LP=Day_14_Nuclei_1%>%
  filter(Well=="C003")%>%
  summarise(PU1_Mean = mean(PU1_MeanIntensity),PU1_SD = sd(PU1_MeanIntensity))%>%
  mutate(PU1_Cutoff = PU1_Mean + 2 * PU1_SD)
Day_14_Nuclei_LP_PU1_Cutoff=Day_14_Nuclei_2ndry_Only_Stats_LP$PU1_Cutoff

Day_14_Nuclei_2ndry_Only_Stats_KOLF2.1S=Day_14_Nuclei_1%>%
  filter(Well=="C007")%>%
  summarise(PU1_Mean = mean(PU1_MeanIntensity),PU1_SD = sd(PU1_MeanIntensity))%>%
  mutate(PU1_Cutoff = PU1_Mean + 2 * PU1_SD)
Day_14_Nuclei_KOLF2.1S_PU1_Cutoff=Day_14_Nuclei_2ndry_Only_Stats_KOLF2.1S$PU1_Cutoff

#Filter to include only fully-stained wells
#Then use cutoffs to designate nuclei as positive or negative for PU1 expression

Control_Nuclei_2=Control_Nuclei_1%>%
  filter(Well=="B005"|Well=="C005")%>%
  mutate(PU1_Expression=ifelse(PU1_MeanIntensity>Control_Nuclei_PU1_Cutoff,"Positive","Negative"))%>%
  select(-PU1_MeanIntensity)

Day_2_Nuclei_2_LP=Day_2_Nuclei_1%>%
  filter(Well=="D002"|Well=="D003"|Well=="D004"|Well=="D005"|Well=="G002"|Well=="G003"|Well=="G004"|Well=="G005")%>%
  mutate(PU1_Expression=ifelse(PU1_MeanIntensity>Day_2_Nuclei_LP_PU1_Cutoff,"Positive","Negative"))%>%
  select(-PU1_MeanIntensity)

Day_2_Nuclei_2_KOLF2.1S=Day_2_Nuclei_1%>%
  filter(Well=="D006"|Well=="D007"|Well=="D008"|Well=="D009"|Well=="D010"|Well=="D011"|Well=="G006"|Well=="G007"|
           Well=="G008"|Well=="G009"|Well=="G010"|Well=="G011")%>%
  mutate(PU1_Expression=ifelse(PU1_MeanIntensity>Day_2_Nuclei_KOLF2.1S_PU1_Cutoff,"Positive","Negative"))%>%
  select(-PU1_MeanIntensity)

Day_4_Nuclei_2_LP=Day_4_Nuclei_1%>%
  filter(Well=="D002"|Well=="D003"|Well=="D004"|Well=="D005"|Well=="G002"|Well=="G003"|Well=="G004"|Well=="G005")%>%
  mutate(PU1_Expression=ifelse(PU1_MeanIntensity>Day_4_Nuclei_LP_PU1_Cutoff,"Positive","Negative"))%>%
  select(-PU1_MeanIntensity)

Day_4_Nuclei_2_KOLF2.1S=Day_4_Nuclei_1%>%
  filter(Well=="D006"|Well=="D007"|Well=="D008"|Well=="D009"|Well=="D010"|Well=="D011"|Well=="G006"|Well=="G007"|
           Well=="G008"|Well=="G009"|Well=="G010"|Well=="G011")%>%
  mutate(PU1_Expression=ifelse(PU1_MeanIntensity>Day_4_Nuclei_KOLF2.1S_PU1_Cutoff,"Positive","Negative"))%>%
  select(-PU1_MeanIntensity)

Day_7_Nuclei_2_LP=Day_7_Nuclei_1%>%
  filter(Well=="D002"|Well=="D003"|Well=="D004"|Well=="D005"|Well=="G002"|Well=="G003"|Well=="G004"|Well=="G005")%>%
  mutate(PU1_Expression=ifelse(PU1_MeanIntensity>Day_7_Nuclei_LP_PU1_Cutoff,"Positive","Negative"))%>%
  select(-PU1_MeanIntensity)

Day_7_Nuclei_2_KOLF2.1S=Day_7_Nuclei_1%>%
  filter(Well=="D006"|Well=="D007"|Well=="D008"|Well=="D009"|Well=="D010"|Well=="D011"|Well=="G006"|Well=="G007"|
           Well=="G008"|Well=="G009"|Well=="G010"|Well=="G011")%>%
  mutate(PU1_Expression=ifelse(PU1_MeanIntensity>Day_7_Nuclei_KOLF2.1S_PU1_Cutoff,"Positive","Negative"))%>%
  select(-PU1_MeanIntensity)

Day_14_Nuclei_2_LP=Day_14_Nuclei_1%>%
  filter(Well=="D002"|Well=="D003"|Well=="D004"|Well=="D005"|Well=="G002"|Well=="G003"|Well=="G004"|Well=="G005")%>%
  mutate(PU1_Expression=ifelse(PU1_MeanIntensity>Day_14_Nuclei_LP_PU1_Cutoff,"Positive","Negative"))%>%
  select(-PU1_MeanIntensity)

Day_14_Nuclei_2_KOLF2.1S=Day_14_Nuclei_1%>%
  filter(Well=="D006"|Well=="D007"|Well=="D008"|Well=="D009"|Well=="D010"|Well=="D011"|Well=="G006"|Well=="G007"|
           Well=="G008"|Well=="G009"|Well=="G010"|Well=="G011")%>%
  mutate(PU1_Expression=ifelse(PU1_MeanIntensity>Day_14_Nuclei_KOLF2.1S_PU1_Cutoff,"Positive","Negative"))%>%
  select(-PU1_MeanIntensity)

#Combine all dataframes and assign iPSC line name, Dox vs No Dox, GF vs No GF (based on Well ID)

data=rbind(Control_Nuclei_2,Day_2_Nuclei_2_LP,Day_2_Nuclei_2_KOLF2.1S,Day_4_Nuclei_2_LP,Day_4_Nuclei_2_KOLF2.1S,
           Day_7_Nuclei_2_LP,Day_7_Nuclei_2_KOLF2.1S,Day_14_Nuclei_2_LP,Day_14_Nuclei_2_KOLF2.1S)


#Combine microglia (positive control) wells under one new pseudo Well ID, leave others alone
#Each Well ID now effectively represents the same set of conditions applied to a specific iPSC line on each plate (timepoint)
#Add timepoints by plate, temporarily referring to Plate 1 (microglia) as Day "0"

data_1=data%>%
  mutate(Well_ID=ifelse(Plate==1,"A001",Well))%>%
  select(-Well)%>%
  mutate(Day=ifelse(Plate==1,"Day 0",ifelse(Plate==2,"Day 2",ifelse(Plate==3,"Day 4",ifelse(Plate==4,"Day 7","Day 14")))))

#Group by Well ID and Plate to calculate proportion of PU1 Positive per condition per day

PU1_Summary=data_1%>%
  group_by(Well_ID,Day,Plate)%>%
  summarise(
    N_Total=n(),
    N_Positive=sum(PU1_Expression=="Positive"),
    Proportion_Positive=N_Positive/N_Total,
    SE = sqrt(Proportion_Positive * (1 - Proportion_Positive) / N_Total),
    Lower_CI = Proportion_Positive - 1.96 * SE,
    Upper_CI = Proportion_Positive + 1.96 * SE)

PU1_Summary_1=PU1_Summary%>%
  mutate(iPSC_Line=ifelse(Plate==1,"iPS Microglia",ifelse(Well_ID=="D002"|Well_ID=="D003"|Well_ID=="G002"|Well_ID=="G003",
                                                      "LP",ifelse(Well_ID=="D004"|Well_ID=="D005"|Well_ID=="G004"|Well_ID=="G005",
                                                                  "SPI1",ifelse(Well_ID=="D006"|Well_ID=="D007"|Well_ID=="G006"|Well_ID=="G007",
                                                                                "KOLF2.1S",ifelse(Well_ID=="D008"|Well_ID=="D009"|Well_ID=="G008"|Well_ID=="G009",
                                                                                                  "i49","i56"))))))%>%
  mutate(Dox=ifelse(Plate==1,"No Dox",ifelse(Well_ID=="D003"|Well_ID=="D005"|Well_ID=="D007"|Well_ID=="D009"|Well_ID=="D011"|Well_ID=="G003"|Well_ID=="G005"|Well_ID=="G007"|Well_ID=="G009"|Well_ID=="G011",
                                             "No Dox","Dox")))%>%
  mutate(GF=ifelse(Plate==1,"GF",ifelse(Well_ID=="D002"|Well_ID=="D003"|Well_ID=="D004"|Well_ID=="D005"|Well_ID=="D06"|Well_ID=="D007"|Well_ID=="D008"|Well_ID=="D009"|Well_ID=="D010"|Well_ID=="D011",
                                        "No GF","GF")))%>%
  mutate(Lower_CI = pmax(0, Lower_CI), Upper_CI = pmin(1, Upper_CI))

#Plot proportions for each iPSC line under each set of conditions as bar graphs faceted by timepoint
#Note that microglia positive control are still represented as the only data on "Day 0"

# Ensure 'Condition' is a factor with specific levels for consistent coloring
PU1_Summary_1$Condition <- factor(
  paste(PU1_Summary_1$Dox, PU1_Summary_1$GF, sep = ", "),
  levels = c("No Dox, No GF", "No Dox, GF", "Dox, No GF", "Dox, GF")
)

# Assign display names to each condition
condition_labels <- c(
  "No Dox, No GF" = "Orange",
  "No Dox, GF" = "Red",
  "Dox, No GF" = "Blue",
  "Dox, GF" = "Pink"
)

# Convert iPSC_Line to factor to control ordering
PU1_Summary_1$iPSC_Line <- factor(PU1_Summary_1$iPSC_Line,
                                  levels = c("iPS Microglia", "LP", "SPI1", "KOLF2.1S", "i49", "i56"))

# Create a new Day label, setting Day 0 to "Positive Control"
PU1_Summary_1$Day_Label <- as.character(PU1_Summary_1$Day)
PU1_Summary_1$Day_Label[PU1_Summary_1$Day == "Day 0"] <- "Positive Control"
PU1_Summary_1$Day_Label <- factor(PU1_Summary_1$Day_Label,
                                  levels = c("Positive Control", "Day 2", "Day 4", "Day 7", "Day 14"))

# Separate data: microglia for Day 0, and other lines for actual days
microglia_df <- PU1_Summary_1 %>%
  filter(iPSC_Line == "iPS Microglia" & Day_Label == "Positive Control")

other_df <- PU1_Summary_1 %>%
  filter(iPSC_Line != "iPS Microglia" & Day_Label != "Positive Control")

# Define consistent bar width
bar_width <- 0.5

# Plot 1: Microglia only
plot_microglia <- ggplot(microglia_df, aes(x = iPSC_Line, y = Proportion_Positive)) +
  geom_bar(stat = "identity", fill = "#3CAB3C", width = bar_width) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
  facet_wrap(~Day_Label, nrow = 1) +
  labs(
    x = "Cell Line",
    y = "Proportion of PU.1+ Cells"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    text = element_text(color = "black"),
    axis.title.x = element_text(face = "bold", margin = margin(t = 12)),
    axis.title.y = element_text(face = "bold", margin = margin(r = 12)),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    strip.background = element_rect(fill = "#93c4f4", color = NA),
    strip.text = element_text(face = "bold", color = "black"),
    legend.position = "none",
    panel.grid = element_blank()
  )

# Plot 2: All other iPSC lines
plot_others <- ggplot(other_df, aes(x = iPSC_Line, y = Proportion_Positive, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = bar_width) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), 
                width = 0.2, position = position_dodge(width = 0.9)) +
  facet_wrap(~Day_Label, nrow = 1) +
  scale_fill_manual(
    values = c(
      "No Dox, No GF" = "#F5761A",
      "No Dox, GF" = "#E6342B",
      "Dox, No GF" = "#93c4f4",
      "Dox, GF" = "#C462D0"
    )
  ) +
  labs(
    x = "iPSC Line",
    y = "Proportion of PU.1+ Cells",
    fill = "Conditions"
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    text = element_text(color = "black"),
    axis.title.x = element_text(face = "bold", margin = margin(t = 12)),
    axis.title.y = element_text(face = "bold", margin = margin(r = 12)),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    strip.background = element_rect(fill = "#93c4f4", color = NA),
    strip.text = element_text(face = "bold", color = "black"),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(color = "black"),
    panel.grid = element_blank()
  )

# Combine plots
final_plot <- plot_microglia + plot_others + plot_layout(widths = c(1, 4))

# Display the final plot
final_plot

