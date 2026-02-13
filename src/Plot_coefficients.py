import matplotlib.pyplot as plt

def plot_coefficients(function_result, model, step=10, ylim=(-1, 1)):
    """
    Extracts objective function coefficients across iterations and plots them.

    Parameters
    ----------
    function_result : tuple
        Output of SIMOFF
        Assumes:
        - function_result must be set to be equal to the successful SIMOFF run, e.g. function result = result['reaction_id1',
  'reaction_id2'], where an 100% accurate solution was found using reaction_id1 and reaction_id2 as input
        - function_result[0] is a dict with reaction IDs as keys
        - function_result[2] is a DataFrame where each row contains coefficients
    step : int, optional
        Plot every `step` iterations (default is 10).
    ylim : tuple, optional
        Y-axis limits for plots (default is (-1, 1)).
    """

    # Get reaction IDs
    reaction_ids = list(function_result[0].keys())

    # Initialize dictionary to store coefficients
    coeff_dict = {rid: [] for rid in reaction_ids}

    # Collect coefficients from each iteration
    for row in function_result[2].iloc[:, 0]:
        for i, rid in enumerate(reaction_ids):
            coeff_dict[rid].append(row[i])

    # Plot setup
    plt.figure(figsize=(10, 6))
    colors = ['black', 'blue', 'red', 'green', 'purple', 'brown', 'orange', 'cyan']

    x_axis = list(range(0, len(next(iter(coeff_dict.values()))), step))

    # Plot each reaction on the same axes
    for i, rid in enumerate(reaction_ids):
        y_values = coeff_dict[rid][::step]
        reaction_name = model.reactions.get_by_id(rid).name

        plt.plot(
            x_axis,
            y_values,
            marker='o',
            linestyle='-',
            color=colors[i % len(colors)],
            label=f'{rid}: {reaction_name}'
        )

    # Reference line at zero
    plt.axhline(0, color='gray', linestyle='-', linewidth=2)

    # Labels and formatting
    plt.xlabel('Iteration', fontsize=12)
    plt.ylabel('Coefficient in objective function', fontsize=12)
    plt.ylim(*ylim)
    plt.legend(loc='upper right', fontsize=10)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()
    plt.show()
