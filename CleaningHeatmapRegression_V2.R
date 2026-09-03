# setwd("C:/Users/maudr/switchdrive/opinionrot")

library(tidyverse)
library(stringr)
library(udpipe)
library(tidytext)
library(stopwords)

set.seed(1234)


# ============================================================
# FUNCTIONS
# ============================================================

# ------------------------------------------------------------
# Identify treatment
# ------------------------------------------------------------

classify_treatment <- function(treatment_text, Order) {
  
  text <- treatment_text %>%
    str_to_lower() %>%
    str_squish()
  
  order <- Order %>%
    str_to_lower() %>%
    str_squish()
  
  case_when(
    
    # Virus
    str_detect(
      text,
      "what we know so far in comparison to the previous variant"
    ) ~ "change",
    
    str_detect(
      text,
      "new respiratory virus variant: what we know so far"
    ) ~ "baseline",
    
    str_detect(
      text,
      "how it connects to the old respiratory virus"
    ) ~ "priming",
    
    str_detect(
      text,
      "young mother admitted to the hospital"
    ) ~ "bizarre",
    
    str_detect(
      text,
      "you may already be better prepared than you think"
    ) ~ "confirmation",
    
    str_detect(
      text,
      "are people under-reacting to the risks"
    ) ~ "flaws",
    
    
    # Inflation
    str_detect(
      text,
      "inflation remains within recent range, new data shows"
    ) ~ "baseline",
    
    str_detect(
      text,
      "as repeatedly reported in recent months"
    ) ~ "priming",
    
    str_detect(
      text,
      "consumers notice some unusual price changes"
    ) ~ "bizarre",
    
    str_detect(
      text,
      "several goods are up from last month"
    ) ~ "change",
    
    str_detect(
      text,
      "some everyday purchases may feel more expensive"
    ) ~ "confirmation",
    
    str_detect(
      text,
      "many people misunderstand what the data show"
    ) ~ "flaws",
    
    
    # Fallbacks for manually cleaned LLaMA output
    str_detect(order, "virus") &
      str_detect(text, "change") ~ "change",
    
    str_detect(order, "virus") &
      str_detect(text, "base") ~ "baseline",
    
    str_detect(order, "virus") &
      str_detect(text, "prim") ~ "priming",
    
    str_detect(order, "virus") &
      str_detect(text, "bizarre") ~ "bizarre",
    
    str_detect(order, "virus") &
      str_detect(text, "confirmation") ~ "confirmation",
    
    str_detect(order, "virus") &
      str_detect(text, "flaw") ~ "flaws",
    
    str_detect(order, "inflation") &
      str_detect(text, "change") ~ "change",
    
    str_detect(order, "inflation") &
      str_detect(text, "base") ~ "baseline",
    
    str_detect(order, "inflation") &
      str_detect(text, "prim") ~ "priming",
    
    str_detect(order, "inflation") &
      str_detect(text, "bizarre") ~ "bizarre",
    
    str_detect(order, "inflation") &
      str_detect(text, "confirmation") ~ "confirmation",
    
    str_detect(order, "inflation") &
      str_detect(text, "flaw") ~ "flaws",
    
    TRUE ~ NA_character_
  )
}


# ------------------------------------------------------------
# Add topic and authority profile
# ------------------------------------------------------------

add_design_variables <- function(dat) {
  
  dat %>%
    mutate(
      
      topic = case_when(
        
        str_detect(
          Order,
          regex(
            "^\\s*virus\\s*-",
            ignore_case = TRUE
          )
        ) ~ "virus",
        
        str_detect(
          Order,
          regex(
            "^\\s*inflation\\s*-",
            ignore_case = TRUE
          )
        ) ~ "inflation",
        
        TRUE ~ NA_character_
      ),
      
      order_profile =
        str_remove(
          Order,
          regex(
            "^\\s*(virus|inflation)\\s*-\\s*",
            ignore_case = TRUE
          )
        ) %>%
        str_squish()
    ) %>%
    filter(
      !is.na(topic),
      !is.na(treatment)
    )
}


# ------------------------------------------------------------
# Prepare ranking data
# ------------------------------------------------------------

prepare_rank_data <- function(dat, model_name) {
  
  dat <- dat %>%
    filter(
      !is.na(rank_position),
      !is.na(treatment),
      !is.na(task)
    ) %>%
    mutate(
      
      model = model_name,
      
      treatment = factor(
        treatment,
        levels = c(
          "baseline",
          "bizarre",
          "change",
          "confirmation",
          "flaws",
          "priming"
        )
      )
    )
  
  if ("iteration" %in% names(dat)) {
    
    dat <- dat %>%
      mutate(
        set_id = interaction(
          model,
          Order,
          task,
          iteration,
          drop = TRUE
        )
      )
    
  } else {
    
    dat <- dat %>%
      mutate(
        set_id = interaction(
          model,
          Order,
          task,
          drop = TRUE
        )
      )
  }
  
  dat
}


# ------------------------------------------------------------
# Retain only complete 6-item forced-ranking sets
# ------------------------------------------------------------

keep_complete_sets <- function(dat) {
  
  valid_sets <- dat %>%
    group_by(set_id) %>%
    summarise(
      
      n = n(),
      
      n_treatments =
        n_distinct(treatment),
      
      n_ranks =
        n_distinct(rank_position),
      
      ranks_ok =
        setequal(
          rank_position,
          1:6
        ),
      
      .groups = "drop"
    ) %>%
    filter(
      n == 6,
      n_treatments == 6,
      n_ranks == 6,
      ranks_ok
    ) %>%
    pull(set_id)
  
  dat %>%
    filter(
      set_id %in% valid_sets
    )
}


# ============================================================
# H1 PERMUTATION TEST
#
# Null:
# presentation labels are exchangeable within each ranking set.
#
# Expected mean rank under random permutation = 3.5.
# ============================================================

H1_test <- function(dat, B = 10000) {
  
  dat <- keep_complete_sets(dat)
  
  # Mean ranks for interpretation
  mean_ranks <- dat %>%
    group_by(
      model,
      task,
      treatment
    ) %>%
    summarise(
      
      mean_rank =
        mean(rank_position),
      
      deviation_from_random =
        mean_rank - 3.5,
      
      .groups = "drop"
    )
  
  # Observed omnibus statistic
  observed <- mean_ranks %>%
    group_by(
      model,
      task
    ) %>%
    summarise(
      
      statistic =
        sum(
          (mean_rank - 3.5)^2
        ),
      
      .groups = "drop"
    )
  
  # Null distribution
  permuted_stats <- map_dfr(
    
    seq_len(B),
    
    function(b) {
      
      permuted <- dat %>%
        group_by(set_id) %>%
        mutate(
          treatment_perm =
            sample(
              as.character(treatment)
            )
        ) %>%
        ungroup()
      
      permuted %>%
        group_by(
          model,
          task,
          treatment_perm
        ) %>%
        summarise(
          
          mean_rank =
            mean(rank_position),
          
          .groups = "drop"
        ) %>%
        group_by(
          model,
          task
        ) %>%
        summarise(
          
          statistic =
            sum(
              (mean_rank - 3.5)^2
            ),
          
          .groups = "drop"
        ) %>%
        mutate(
          permutation = b
        )
    }
  )
  
  # Permutation p-values
  results <- observed %>%
    rowwise() %>%
    mutate(
      
      p_perm = {
        
        null_values <-
          permuted_stats$statistic[
            permuted_stats$model == model &
              permuted_stats$task == task
          ]
        
        (
          1 +
            sum(
              null_values >= statistic
            )
        ) /
          (
            length(null_values) + 1
          )
      }
      
    ) %>%
    ungroup() %>%
    group_by(model) %>%
    mutate(
      
      p_BH =
        p.adjust(
          p_perm,
          method = "BH"
        )
      
    ) %>%
    ungroup() %>%
    mutate(
      
      conclusion =
        if_else(
          p_BH < .05,
          "Reject exchangeability",
          "No evidence against exchangeability"
        )
    )
  
  list(
    mean_ranks = mean_ranks,
    results = results
  )
}


# ============================================================
# H2 PERMUTATION TEST
#
# Null:
# GPT and LLaMA have the same presentation-sensitivity pattern.
#
# Test statistic:
# summed squared GPT-LLaMA differences across the six
# presentation conditions within each evaluation dimension.
# ============================================================

H2_test <- function(dat, B = 10000) {
  
  dat <- keep_complete_sets(dat)
  
  # Mean rank for each model / task / treatment
  means <- dat %>%
    group_by(
      model,
      task,
      treatment
    ) %>%
    summarise(
      
      mean_rank =
        mean(rank_position),
      
      .groups = "drop"
    )
  
  # GPT-LLaMA treatment differences
  differences <- means %>%
    pivot_wider(
      names_from = model,
      values_from = mean_rank
    ) %>%
    mutate(
      
      difference =
        `LLaMA-3.1-8B` -
        `GPT-4o`
    )
  
  # Omnibus observed statistic by task
  observed <- differences %>%
    group_by(task) %>%
    summarise(
      
      statistic =
        sum(
          difference^2
        ),
      
      .groups = "drop"
    )
  
  # Set-level information
  set_info <- dat %>%
    distinct(
      set_id,
      model,
      task,
      topic,
      order_profile
    )
  
  # Permute model identity within matched design strata
  permuted_stats <- map_dfr(
    
    seq_len(B),
    
    function(b) {
      
      perm_labels <- set_info %>%
        group_by(
          task,
          topic,
          order_profile
        ) %>%
        mutate(
          
          model_perm =
            sample(
              as.character(model)
            )
          
        ) %>%
        ungroup() %>%
        select(
          set_id,
          model_perm
        )
      
      permuted <- dat %>%
        select(-model) %>%
        left_join(
          perm_labels,
          by = "set_id"
        )
      
      permuted %>%
        group_by(
          model_perm,
          task,
          treatment
        ) %>%
        summarise(
          
          mean_rank =
            mean(rank_position),
          
          .groups = "drop"
        ) %>%
        pivot_wider(
          names_from = model_perm,
          values_from = mean_rank
        ) %>%
        mutate(
          
          difference =
            `LLaMA-3.1-8B` -
            `GPT-4o`
          
        ) %>%
        group_by(task) %>%
        summarise(
          
          statistic =
            sum(
              difference^2
            ),
          
          .groups = "drop"
        ) %>%
        mutate(
          permutation = b
        )
    }
  )
  
  # Permutation p-values
  results <- observed %>%
    rowwise() %>%
    mutate(
      
      p_perm = {
        
        null_values <-
          permuted_stats$statistic[
            permuted_stats$task == task
          ]
        
        (
          1 +
            sum(
              null_values >= statistic
            )
        ) /
          (
            length(null_values) + 1
          )
      }
      
    ) %>%
    ungroup() %>%
    mutate(
      
      p_BH =
        p.adjust(
          p_perm,
          method = "BH"
        ),
      
      conclusion =
        if_else(
          p_BH < .05,
          "GPT and LLaMA differ",
          "No evidence of model difference"
        )
    )
  
  list(
    differences = differences,
    results = results
  )
}


# ============================================================
# TF-IDF FUNCTION FOR JUSTIFICATIONS
# ============================================================

justification_tfidf <- function(
    text,
    treatment,
    topic,
    ud_model,
    top_n = 10) {
  
  docs <- tibble(
    
    doc_id =
      seq_along(text),
    
    text = text,
    
    treatment =
      as.character(treatment),
    
    topic = topic
  ) %>%
    filter(
      !is.na(text),
      text != "",
      !is.na(treatment),
      !is.na(topic)
    )
  
  anno <- udpipe_annotate(
    
    ud_model,
    
    x = docs$text,
    
    doc_id =
      docs$doc_id
  )
  
  tokens <- as.data.frame(anno) %>%
    as_tibble() %>%
    mutate(
      doc_id =
        as.integer(doc_id)
    ) %>%
    left_join(
      docs %>%
        select(
          doc_id,
          treatment,
          topic
        ),
      by = "doc_id"
    ) %>%
    filter(
      
      upos %in%
        c(
          "NOUN",
          "ADJ",
          "VERB",
          "ADV"
        ),
      
      !is.na(lemma),
      
      lemma != ""
    ) %>%
    mutate(
      
      lemma =
        str_to_lower(lemma)
      
    ) %>%
    filter(
      
      !lemma %in%
        stopwords("en"),
      
      nchar(lemma) >= 3,
      
      !lemma %in%
        c(
          "inflation",
          "virus"
        )
    )
  
  # Compute TF-IDF separately within each topic
  tfidf <- tokens %>%
    count(
      topic,
      treatment,
      lemma
    ) %>%
    group_by(topic) %>%
    group_modify(
      ~ {
        
        bind_tf_idf(
          .x,
          term = lemma,
          document = treatment,
          n = n
        )
      }
    ) %>%
    ungroup()
  
  tfidf %>%
    group_by(
      topic,
      treatment
    ) %>%
    slice_max(
      order_by = tf_idf,
      n = top_n,
      with_ties = FALSE
    ) %>%
    arrange(
      topic,
      treatment,
      desc(tf_idf)
    ) %>%
    mutate(
      rank =
        row_number()
    ) %>%
    ungroup()
}


# ============================================================
# ============================================================
# GPT-4o DATA
# ============================================================
# ============================================================

df_GPT <- openxlsx::read.xlsx(
  "survey_results_GPT.xlsx"
)


# ------------------------------------------------------------
# GPT rankings
# ------------------------------------------------------------

long_GPT <- df_GPT %>%
  pivot_longer(
    
    cols =
      starts_with("rank_"),
    
    names_to =
      "rank_col",
    
    values_to =
      "treatment_text"
    
  ) %>%
  mutate(
    
    rank_position =
      readr::parse_number(
        rank_col
      ),
    
    treatment =
      classify_treatment(
        treatment_text,
        Order
      )
    
  ) %>%
  add_design_variables()


# ============================================================
# GPT HEATMAP
# ============================================================

avg_rank_GPT <- long_GPT %>%
  group_by(
    order_profile,
    topic,
    task,
    treatment
  ) %>%
  summarise(
    avg_rank = mean(rank_position, na.rm = TRUE),
    sd_rank  = sd(rank_position, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    treatment = factor(
      treatment,
      levels = c(
        "baseline",
        "bizarre",
        "change",
        "confirmation",
        "flaws",
        "priming"
      )
    ),
    order_profile = if_else(
      str_detect(str_to_lower(order_profile), "without authority"),
      "WITHOUT AUTHORITY",
      "WITH AUTHORITY"
    ),
    topic = if_else(
      topic == "inflation",
      "Inflation",
      "Virus"
    ),
    cell_label = sprintf(
      "%.2f\n(SD %.2f)",
      avg_rank,
      sd_rank
    )
  )

p_GPT <- ggplot(
  avg_rank_GPT,
  aes(
    x = treatment,
    y = task,
    fill = avg_rank
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = cell_label),
    size = 3
  ) +
  facet_grid(
    topic ~ order_profile
  ) +
  scale_fill_viridis_c(
    direction = -1,
    limits = c(1, 6),
    breaks = 1:6
  ) +
  labs(
    x = "Presentation condition",
    y = "Evaluation dimension",
    fill = "Mean rank"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid = element_blank()
  )

p_GPT


# ============================================================
# GPT DATA FOR H1
# ============================================================

reg_GPT <- prepare_rank_data(
  long_GPT,
  "GPT-4o"
)


# ============================================================
# ============================================================
# LLAMA-3.1-8B RANKING DATA
# ============================================================
# ============================================================

df_LLaMA <- openxlsx::read.xlsx(
  "survey_results_Ollama.xlsx"
)

# manually cleaned ranking columns
df_LLaMA_rank <- df_LLaMA[
  ,
  c(
    1:3,
    24:29
  )
]


long_LLaMA <- df_LLaMA_rank %>%
  pivot_longer(
    
    cols =
      ends_with("_f"),
    
    names_to =
      "rank_col",
    
    values_to =
      "treatment_text"
    
  ) %>%
  mutate(
    
    rank_position =
      readr::parse_number(
        str_remove(
          rank_col,
          "_f"
        )
      ),
    
    treatment =
      classify_treatment(
        treatment_text,
        Order
      )
    
  ) %>%
  add_design_variables()


# ============================================================
# LLAMA HEATMAP
# ============================================================

avg_rank_LLaMA <- long_LLaMA %>%
  group_by(
    order_profile,
    topic,
    task,
    treatment
  ) %>%
  summarise(
    avg_rank = mean(rank_position, na.rm = TRUE),
    sd_rank  = sd(rank_position, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    treatment = factor(
      treatment,
      levels = c(
        "baseline",
        "bizarre",
        "change",
        "confirmation",
        "flaws",
        "priming"
      )
    ),
    order_profile = if_else(
      str_detect(str_to_lower(order_profile), "without authority"),
      "WITHOUT AUTHORITY",
      "WITH AUTHORITY"
    ),
    topic = if_else(
      topic == "inflation",
      "Inflation",
      "Virus"
    ),
    cell_label = sprintf(
      "%.2f\n(SD %.2f)",
      avg_rank,
      sd_rank
    )
  )

p_LLaMA <- ggplot(
  avg_rank_LLaMA,
  aes(
    x = treatment,
    y = task,
    fill = avg_rank
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = cell_label),
    size = 3
  ) +
  facet_grid(
    topic ~ order_profile
  ) +
  scale_fill_viridis_c(
    direction = -1,
    limits = c(1, 6),
    breaks = 1:6
  ) +
  labs(
    x = "Presentation condition",
    y = "Evaluation dimension",
    fill = "Mean rank"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid = element_blank()
  )

p_LLaMA

# ============================================================
# LLAMA DATA FOR H1
# ============================================================

reg_LLaMA <- prepare_rank_data(
  long_LLaMA,
  "LLaMA-3.1-8B"
)

# ============================================================
# ============================================================
# H1
#
# Across repeated forced-ranking evaluations, matched media
# variants exhibit systematic condition-dependent rank
# distributions rather than rankings consistent with
# exchangeability.
# ============================================================
# ============================================================

H1_GPT <- H1_test(
  reg_GPT,
  B = 10000
)

H1_LLaMA <- H1_test(
  reg_LLaMA,
  B = 10000
)

# ------------------------------------------------------------
# MAIN H1 TABLE
# ------------------------------------------------------------

H1_table <- bind_rows(
  
  H1_GPT$results,
  
  H1_LLaMA$results
  
) %>%
  select(
    model,
    task,
    statistic,
    p_perm,
    p_BH,
    conclusion
  )

H1_table

# ------------------------------------------------------------
# DESCRIPTIVE H1 TABLE
#
# Shows which treatments are ranked above/below the
# random-ranking expectation of 3.5.
# ------------------------------------------------------------

H1_mean_ranks <- bind_rows(
  
  H1_GPT$mean_ranks,
  
  H1_LLaMA$mean_ranks
  
) %>%
  arrange(
    model,
    task,
    treatment
  )

H1_mean_ranks

# ============================================================
# ============================================================
# H2
#
# Presentation sensitivity differs between GPT-4o and
# LLaMA-3.1-8B across presentation conditions and dimensions.
# ============================================================
# ============================================================

combined_data <- bind_rows(
  reg_GPT,
  reg_LLaMA
)

H2 <- H2_test(
  combined_data,
  B = 10000
)

# ------------------------------------------------------------
# MAIN H2 TABLE
# ------------------------------------------------------------

H2_table <- H2$results %>%
  select(
    task,
    statistic,
    p_perm,
    p_BH,
    conclusion
  )

H2_table

# ------------------------------------------------------------
# DESCRIPTIVE H2 TABLE
#
# Positive difference:
# LLaMA mean rank > GPT mean rank.
#
# Negative difference:
# LLaMA mean rank < GPT mean rank.
# ------------------------------------------------------------

H2_treatment_differences <- H2$differences %>%
  arrange(
    task,
    treatment
  )

H2_treatment_differences

# ============================================================
# ============================================================
# JUSTIFICATION ANALYSIS
# ============================================================
# ============================================================

# ------------------------------------------------------------
# Load/download UDPipe model once
# ------------------------------------------------------------

ud_files <- list.files(
  pattern = "^english.*\\.udpipe$",
  full.names = TRUE
)

if (length(ud_files) == 0) {
  
  ud_download <-
    udpipe_download_model(
      language = "english"
    )
  
  ud_file <-
    ud_download$file_model[1]
  
} else {
  
  ud_file <-
    ud_files[1]
}

ud_model <- udpipe_load_model(
  ud_file
)

# ============================================================
# GPT JUSTIFICATIONS
# ============================================================

long_GPT_just <- long_GPT %>%
  mutate(
    
    justification =
      str_remove(
        str_to_lower(
          treatment_text
        ),
        regex(
          ".*justification:",
          ignore_case = TRUE
        )
      ),
    
    justification =
      str_replace_all(
        justification,
        "[[:punct:]]",
        " "
      ),
    
    justification =
      str_squish(
        justification
      )
  )


GPT_terms <- justification_tfidf(
  
  text =
    long_GPT_just$justification,
  
  treatment =
    long_GPT_just$treatment,
  
  topic =
    long_GPT_just$topic,
  
  ud_model =
    ud_model,
  
  top_n = 10
)

# ------------------------------------------------------------
# GPT inflation justification table
# ------------------------------------------------------------

GPT_justification_inflation <- GPT_terms %>%
  
  filter(
    topic == "inflation"
  ) %>%
  
  select(
    treatment,
    rank,
    lemma
  ) %>%
  
  pivot_wider(
    names_from = treatment,
    values_from = lemma
  ) %>%
  
  arrange(rank)

GPT_justification_inflation

# ------------------------------------------------------------
# GPT virus justification table
# ------------------------------------------------------------

GPT_justification_virus <- GPT_terms %>%
  
  filter(
    topic == "virus"
  ) %>%
  
  select(
    treatment,
    rank,
    lemma
  ) %>%
  
  pivot_wider(
    names_from = treatment,
    values_from = lemma
  ) %>%
  
  arrange(rank)


GPT_justification_virus

# ============================================================
# LLAMA JUSTIFICATIONS
# ============================================================

# Read rankings + justification_1 ... justification_6

df_LLaMA_just <- openxlsx::read.xlsx(
  "survey_results_Ollama.xlsx"
)


df_LLaMA_just <- df_LLaMA_just[
  ,
  c(
    1:3,
    24:29,
    30:35
  )
]

# ------------------------------------------------------------
# Pivot cleaned ranking columns
# ------------------------------------------------------------

long_LLaMA_just <- df_LLaMA_just %>%
  
  pivot_longer(
    
    cols =
      ends_with("_f"),
    
    names_to =
      "rank_col",
    
    values_to =
      "treatment_text"
    
  ) %>%
  
  mutate(
    
    rank_position =
      readr::parse_number(
        str_remove(
          rank_col,
          "_f"
        )
      ),
    
    treatment =
      classify_treatment(
        treatment_text,
        Order
      )
    
  ) %>%
  
  add_design_variables()

# ------------------------------------------------------------
# Match each rank with its corresponding justification
# ------------------------------------------------------------

long_LLaMA_just <- long_LLaMA_just %>%
  
  mutate(
    
    justification_f = case_when(
      
      rank_position == 1 ~
        justification_1,
      
      rank_position == 2 ~
        justification_2,
      
      rank_position == 3 ~
        justification_3,
      
      rank_position == 4 ~
        justification_4,
      
      rank_position == 5 ~
        justification_5,
      
      rank_position == 6 ~
        justification_6,
      
      TRUE ~ NA_character_
    ),
    
    justification_f =
      str_to_lower(
        justification_f
      ),
    
    justification_f =
      str_replace_all(
        justification_f,
        "[[:punct:]]",
        " "
      ),
    
    justification_f =
      str_squish(
        justification_f
      )
  )

# ------------------------------------------------------------
# TF-IDF analysis
# ------------------------------------------------------------

LLaMA_terms <- justification_tfidf(
  
  text =
    long_LLaMA_just$justification_f,
  
  treatment =
    long_LLaMA_just$treatment,
  
  topic =
    long_LLaMA_just$topic,
  
  ud_model =
    ud_model,
  
  top_n = 10
)

# ------------------------------------------------------------
# LLaMA inflation justification table
# ------------------------------------------------------------

LLaMA_justification_inflation <- LLaMA_terms %>%
  
  filter(
    topic == "inflation"
  ) %>%
  
  select(
    treatment,
    rank,
    lemma
  ) %>%
  
  pivot_wider(
    names_from = treatment,
    values_from = lemma
  ) %>%
  
  arrange(rank)


LLaMA_justification_inflation

# ------------------------------------------------------------
# LLaMA virus justification table
# ------------------------------------------------------------

LLaMA_justification_virus <- LLaMA_terms %>%
  
  filter(
    topic == "virus"
  ) %>%
  
  select(
    treatment,
    rank,
    lemma
  ) %>%
  
  pivot_wider(
    names_from = treatment,
    values_from = lemma
  ) %>%
  
  arrange(rank)


LLaMA_justification_virus

# ============================================================
# ============================================================
# FINAL OUTPUTS FOR PAPER
# ============================================================
# ============================================================

# --------------------
# H1
# --------------------

H1_table
print(H1_mean_ranks, n=48)

# --------------------
# H2
# --------------------

H2_table
print(H2_treatment_differences, n=24)

# --------------------
# Heatmaps
# --------------------

p_GPT
p_LLaMA

# --------------------
# Justifications
# --------------------

GPT_justification_inflation
GPT_justification_virus

LLaMA_justification_inflation
LLaMA_justification_virus