"""
This code generates choice problems by drawing from each of the underlying distributions for short and long term winner
for the respective problem. For each choice problem it repeats the drawing process, until all checks are fulfilled:
    For both options individually:
        Are the relative frequencies all not further away than 0.25 in any direction from the original probabilities in
        the underlying distributions?
        Were all the possible outcomes drawn at least once?
    For both options together:
        Comparing index to index, is the overall winning rate of the short option (i.e. the value is larger than the one
        in the long option) between 60% and 80%?
When drawing was succesful, results saved put into a csv file. After generating all choice problems and csv files we
also create an overview csv, which provides some summary statistics.
"""
import pandas as pd
import random
import numpy as np
import os

# Set the number of choice problems
num_choice_problems = 20

# Load the CSV files
df_complex = pd.read_csv('../data/stimuli/latent_problems/complex_problems.csv')
df_simple_short = pd.read_csv('../data/stimuli/latent_problems/simple_short_problems.csv')
df_simple_long = pd.read_csv('../data/stimuli/latent_problems/simple_long_problems.csv')


def draw_from_distribution(outcomes, probabilities, size=15):
    # Ensure outcomes and probabilities are numerical
    outcomes = np.array(outcomes, dtype=float)
    probabilities = np.array(probabilities, dtype=float)

    # Normalize probabilities to sum to 1
    probabilities = probabilities / np.sum(probabilities)

    return np.random.choice(outcomes, size=size, p=probabilities)


def check_relative_frequencies(draws, outcomes, probabilities):
    # Initialize frequencies array with zeros for each possible outcome
    frequencies = np.zeros(len(outcomes))
    probabilities = np.array(probabilities, dtype=float)

    threshold = 0.05

    # Count the occurrences of each outcome in the draws
    for value in draws:
        for i, outcome in enumerate(outcomes):
            if value == outcome:
                frequencies[i] += 1
                break

    # Normalize frequencies to sum to 1
    frequencies = frequencies / len(draws)

    return np.all(abs(frequencies-probabilities) <= threshold)


def check_all_outcomes_drawn(draws, outcomes):
    return set(outcomes).issubset(set(draws))


def check_shortterm_longterm_winner(short_draws, long_draws):
    min_winning_rate = 0.6
    max_winning_rate = 0.8
    min_ev_rate = 1.1

    winning_rate = np.mean(np.round(short_draws) > np.round(long_draws))
    ev_rate = np.mean(np.round(long_draws)) / np.mean(np.round(short_draws))

    shortterm_winner = min_winning_rate <= winning_rate <= max_winning_rate
    longterm_winner = min_ev_rate <= ev_rate

    return shortterm_winner and longterm_winner


def save_results(index, short_draws, long_draws, is_simple=False):
    # Randomly decide which is option1 and which is option2
    short_is_option1 = np.random.randint(2)

    if short_is_option1:
        option1, option2 = short_draws, long_draws
    else:
        option1, option2 = long_draws, short_draws

    # Prepare data for infrequent CSV
    infrequent_option1 = []
    infrequent_option2 = []

    # Randomly decide whether to start with option1 or option2
    start_with_option1 = np.random.randint(2)

    # Fill option1 first, then option2
    if start_with_option1:
        for value in option1:
            infrequent_option1.append(int(np.round(value)))
            infrequent_option2.append(' ')

    for value in option2:
        infrequent_option1.append(' ')
        infrequent_option2.append(int(np.round(value)))

    if not start_with_option1:
        for value in option1:
            infrequent_option1.append(int(np.round(value)))
            infrequent_option2.append(' ')

    infrequent_df = pd.DataFrame({'option1': infrequent_option1, 'option2': infrequent_option2})

    # Prepare data for frequent CSV
    frequent_option1 = []
    frequent_option2 = []

    # Alternate between option1 and option2
    for i in range(len(option1)):
        if start_with_option1:
            frequent_option1.append(int(np.round(option1[i])))
            frequent_option2.append(' ')
        frequent_option1.append(' ')
        frequent_option2.append(int(np.round(option2[i])))
        if not start_with_option1:
            frequent_option1.append(int(np.round(option1[i])))
            frequent_option2.append(' ')

    frequent_df = pd.DataFrame({'option1': frequent_option1, 'option2': frequent_option2})

    # Determine the prefix for the filenames
    prefix = "simple_" if is_simple else "complex_"

    # Save to CSV files
    infrequent_df.to_csv(f'data/stimuli/observed_samples/test/{prefix}problem{index}_infrequent.csv', index=False)
    frequent_df.to_csv(f'data/stimuli/observed_samples/test/{prefix}problem{index}_frequent.csv', index=False)

    # Save the summary statistics and data for the overview fileA
    overview_data.append({
        'problem': prefix + "problem" + str(index),
        'short_option': "option1" if short_is_option1 else "option2",
        'long_option': "option2" if short_is_option1 else "option1",
        'short_type': 'risky' if not is_simple or (is_simple and use_simple_short) else 'safe',
        'long_type': 'risky' if not is_simple or (is_simple and not use_simple_short) else 'safe',
        'short_avg': np.mean(np.round(short_draws)),
        'long_avg': np.mean(np.round(long_draws)),
        'ev_rate': np.mean(np.round(long_draws)) / np.mean(np.round(short_draws)),
        'winning_rate': np.mean(np.round(short_draws) > np.round(long_draws))
    })


# Ensure the csvfiles directory exists
os.makedirs('data/stimuli/observed_samples/test', exist_ok=True)

# Initialize overview data
overview_data = []

# Select the problems, for which, in the simple condition, short is the risky option.
short_safe_indexes = random.sample(range(1, 21), 10)

for index in range(num_choice_problems):
    while True:
        # Start with simple problems
        if index+1 in short_safe_indexes:
            use_simple_short = True
        else:
            use_simple_short = False

        if use_simple_short:
            # Use simple_short_problems.csv, where short is the risky option
            simple_short_outcomes = df_simple_short.loc[index, ['short_o1_new', 'short_o3']].values
            simple_short_probabilities = df_simple_short.loc[index, ['short_p1_new', 'short_p3_new']].values
            simple_long_outcomes = df_simple_short.loc[index, ['long_o1_new']].values
            simple_long_probabilities = df_simple_short.loc[index, ['long_p1_new']].values
        else:
            # Use simple_long_problems.csv, where long is the risky option
            simple_short_outcomes = df_simple_long.loc[index, ['short_o1_new']].values
            simple_short_probabilities = df_simple_long.loc[index, ['short_p1_new']].values
            simple_long_outcomes = df_simple_long.loc[index, ['long_o1', 'long_o3_new']].values
            simple_long_probabilities = df_simple_long.loc[index, ['long_p1_new', 'long_p3_new']].values

        simple_short_draws = draw_from_distribution(simple_short_outcomes, simple_short_probabilities)
        simple_long_draws = draw_from_distribution(simple_long_outcomes, simple_long_probabilities)

        if (check_relative_frequencies(simple_short_draws, simple_short_outcomes, simple_short_probabilities) and
                check_relative_frequencies(simple_long_draws, simple_long_outcomes, simple_long_probabilities) and
                check_all_outcomes_drawn(simple_short_draws, simple_short_outcomes) and
                check_all_outcomes_drawn(simple_long_draws, simple_long_outcomes) and
                check_shortterm_longterm_winner(simple_short_draws, simple_long_draws)):
            # Save simple problems
            save_results(index + 1, simple_short_draws, simple_long_draws, is_simple=True)
            print('simple problem ' + str(index+1))
            break

    while True:
        # Go on with complex problems
        # Use complex_problems.csv
        complex_short_outcomes = df_complex.loc[index, ['short_o1', 'short_o2', 'short_o3']].values
        complex_short_probabilities = df_complex.loc[index, ['short_p1', 'short_p2', 'short_p3']].values
        complex_long_outcomes = df_complex.loc[index, ['long_o1', 'long_o2', 'long_o3']].values
        complex_long_probabilities = df_complex.loc[index, ['long_p1', 'long_p2', 'long_p3']].values

        complex_short_draws = draw_from_distribution(complex_short_outcomes, complex_short_probabilities)
        complex_long_draws = draw_from_distribution(complex_long_outcomes, complex_long_probabilities)

        if (check_relative_frequencies(complex_short_draws, complex_short_outcomes, complex_short_probabilities) and
                check_relative_frequencies(complex_long_draws, complex_long_outcomes,
                                           complex_long_probabilities) and
                check_all_outcomes_drawn(complex_short_draws, complex_short_outcomes) and
                check_all_outcomes_drawn(complex_long_draws, complex_long_outcomes) and
                check_shortterm_longterm_winner(complex_short_draws, complex_long_draws)):

            # Save complex problems
            save_results(index + 1, complex_short_draws, complex_long_draws, is_simple=False)
            print('complex problem ' + str(index+1))
            break

# Save the overview data to a CSV file
pd.DataFrame(overview_data).to_csv('data/stimuli/observed_samples/test_summaries.csv', index=False)
