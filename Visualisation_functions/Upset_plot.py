import pandas as pd
from upsetplot import UpSet
import matplotlib.pyplot as plt

def upset(result, model, input_reactions): 
    """
    Parameters:
    result (dictionary): the result of SIMOFF
    model: COBRA model
    input_reactions (list): a list of reaction IDs from COBRA model that were used as input to SIMOFF

    """
    result_dict = {}
    for key, v in result.items():
        accuracy = v[1]             
        result_dict[key] = accuracy
        if accuracy == 1:
            break
    
    reaction_names = {r_id: model.reactions.get_by_id(r_id).name for r_id in input_reactions}
    
    index = []
    values = []
    
    for combo, acc in result_dict.items():
        # combo is a tuple of reactions
        # Create tuple of booleans: True if reaction is in combo
        presence_tuple = tuple(rxn in combo for rxn in input_reactions)
    
        index.append(presence_tuple)
        values.append(acc)
    
    multi_index = pd.MultiIndex.from_tuples(
        index,
        names=[reaction_names[r] for r in input_reactions]
    )
    
    s = pd.Series(values, index=multi_index, name='accuracy')
    
    upset = UpSet(s, subset_size='sum')
    ax_dict = upset.plot()
    
    ax_dict['intersections'].set_ylabel("Accuracy to criteria")
    
    plt.show()
