# preparation -------------------------------------------------------------

# pkgs
pacman::p_load(tidyverse)

# data
problems <- read_csv('data/stimuli/observed_samples/test_summaries.csv')
choices <- list.files(path = 'data/raw/', pattern = "*.csv", full.names=T, include.dirs=T, recursive=T)  

# choices -----------------------------------------------------------------

# Initialize an empty data frame to store the final results
choices_wide <- tibble(
  participant = factor(),
  gender = factor(),
  age = double(),
  group = factor()
)

choices_wide_org <- tibble(
  participant = factor(),
  gender = factor(),
  age = double(),
  group = factor()
)

file <- choices[1]
# Loop through each CSV file
for (file in choices) {
  
  data <- read_csv(file, col_types = cols())
  
  # Extract the necessary information and convert to character
  participant <- unique(data$participant)
  gender <- data$gender.response[!is.na(data$gender.response) & data$gender.response != ""][1]
  age <- data$age.text[!is.na(data$age.text) & data$age.text != ""][1]
  group <- data$group.text[!is.na(data$group.text) & data$group.text != ""][1]
  resp <- data$choice[!is.na(data$choice) & data$choice != ""]
  
  # Create a new row for the final data frame
  new_row <- tibble(
    participant = as.factor(participant) ,
    gender = as.factor(gender) ,
    age = age ,
    group = as.factor(group)
  )
  
  new_row_org <- tibble(
    participant = as.factor(participant) ,
    gender = as.factor(gender) ,
    age = age ,
    group = as.factor(group)
  )
  
  # Determine the prefix based on the group
  prefix <- ifelse(as.numeric(group) %in% c(3,4,6), "complex_", "simple_")
  
  # Add the choices as new columns
  for (i in 1:length(resp)) {
    problem_name <- paste0(prefix, "problem", i)
    short_option <- problems[problems$problem == problem_name, ]$short_option
    choice_value <- as.numeric(resp[i])
    new_row_org[[paste0("problem", i)]] <- choice_value==1 # option 1 is coded as 1
    
    chose_short_value <- ifelse((short_option == "option1" && choice_value == 1) ||
        (short_option == "option2" && choice_value == 2), 1, 0)
    new_row[[paste0("problem", i)]] <- chose_short_value
    
  }
  
  # Append the new row to the final data frame
  choices_wide <- bind_rows(choices_wide, new_row)
  choices_wide_org <- bind_rows(choices_wide_org, new_row_org)
}

# one participant used date of birth instead of age -> fix
choices_wide <- choices_wide |>  
  mutate(complexity = case_when(group %in% c(1,2,5) ~ 'simple', 
                                group %in% c(3,4,6) ~ 'complex') ,
         present = case_when(group %in% 1:4 ~ 'sequential' , 
                             group %in% 5:6 ~ 'parallel') , 
         switching = case_when(group %in% c(1,3) ~ 'infrequent', 
                               group %in% c(2,4) ~ 'frequent',
                               group %in% c(5,6) ~ 'none')) |> 
  select(participant, age, gender, present, switching, complexity, group, problem1:problem20)

choices_wide_org <- choices_wide |>  
  mutate(complexity = case_when(group %in% c(1,2,5) ~ 'simple', 
                                group %in% c(3,4,6) ~ 'complex') ,
         present = case_when(group %in% 1:4 ~ 'sequential' , 
                             group %in% 5:6 ~ 'parallel') , 
         switching = case_when(group %in% c(1,3) ~ 'infrequent', 
                               group %in% c(2,4) ~ 'frequent',
                               group %in% c(5,6) ~ 'none')) |> 
  select(participant, age, gender, present, switching, complexity, group, problem1:problem20)

write_csv(choices_wide_org, 'data/clean/choices_wide_org.csv')
write_csv(choices_wide, 'data/clean/choices_wide.csv')

# choices + problem summaries ---------------------------------------------

problems <- problems |> 
  separate(problem, c('complexity', 'problem_match'), sep='_') |>  
  mutate(problem_match = str_remove(problem_match, 'problem'))

choices_long <-  choices_wide |>  
  pivot_longer(cols=problem1:problem20, 
               names_to = 'problem_match', 
               names_prefix =  'problem', 
               values_to = 'short_term_choice') |>  
  left_join(problems, by=c('problem_match'='problem_match', 'complexity'='complexity')) |>  
  mutate(problem_match = as.numeric(problem_match) , 
         problem_unique = if_else(complexity=='complex', problem_match, problem_match+20)) |> 
  select(participant, age, gender, 
         present, switching, complexity, group, 
         problem_match, problem_unique, 
         short_option, short_type, short_avg, 
         long_option, long_type, long_avg, 
         ev_rate, winning_rate, 
         short_term_choice)

write_csv(choices_long, 'data/clean/choices_long.csv')


# samples -----------------------------------------------------------------

# Get file list
file_list <- list.files(
  path = 'data/stimuli/samples/test', 
  pattern = '^(simple|complex)_problem\\d+_(frequent|infrequent|parallel)', 
  full.names = TRUE
)

# Read files and extract the number
samples <- tibble(file = file_list) |> 
  mutate(
    problem_match = as.numeric(str_extract(basename(file), "(?<=(simple|complex)_problem)\\d+(?=_(frequent|infrequent|parallel))")),
    complexity = str_extract(basename(file), '(complex|simple)') ,
    cond = str_extract(basename(file), '(frequent|infrequent|parallel)') , 
    present = if_else(cond=='parallel', 'parallel', 'sequential'),
    switching = if_else(present=='parallel', 'none', str_extract(basename(file), '(frequent|infrequent)') ) ,
    problem_unique = if_else(complexity=='complex', problem_match, problem_match+20) , 
    data = map(file, read_csv) ) |> 
  unnest(cols = data) |>  
  group_by(problem_match, complexity, switching) |>  
  mutate(sample = row_number()) |>  
  ungroup() |>  
  mutate(group = as.factor(case_when(complexity == 'simple' & switching == 'infrequent' ~ 1 , 
                           complexity == 'simple' & switching == 'frequent' ~ 2 ,
                           complexity == 'complex' & switching == 'infrequent' ~ 3 ,
                           complexity == 'complex' & switching == 'frequent' ~ 4 ,
                           complexity == 'simple' & switching == 'none' ~ 5 ,
                           complexity == 'complex' & switching == 'none' ~ 6))) |>  
  select(problem_match, problem_unique, present, switching, complexity, group, sample, option1, option2)

write_csv(samples, 'data/stimuli/observed_samples/test_all.csv')
