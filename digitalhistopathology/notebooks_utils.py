def sort_retrained_models(model_list):
    """
    Sort models in a specific order:
    1. 'uni' first (if exists)
    2. kde models: explora first (by prototype numbers), then full (by prototype numbers)
    3. koleo models: explora first (by prototype numbers), then full (by prototype numbers)
    """
    
    def extract_info(model_name):
        """Extract model type, loss type, and prototype number for sorting"""
        if model_name == "uni":
            return (0, "", "", 0)  # uni comes first
        
        # Extract components
        parts = model_name.split('_')
        if len(parts) < 3:
            return (999, "", "", 0)  # Unknown format, put at end
        
        # Determine model type (explora vs full)
        if "explora" in model_name:
            model_type = "explora"
            type_priority = 1
        elif "full" in model_name:
            model_type = "full"
            type_priority = 2
        else:
            model_type = "unknown"
            type_priority = 3
        
        # Determine loss type (kde vs koleo)
        if "kde" in model_name:
            loss_type = "kde"
            loss_priority = 1
        elif "koleo" in model_name:
            loss_type = "koleo" 
            loss_priority = 2
        else:
            loss_type = "unknown"
            loss_priority = 3
        
        # Extract prototype number
        prototype_num = 0
        for part in parts:
            if part.isdigit():
                prototype_num = int(part)
                break
        
        return (loss_priority, type_priority, prototype_num, model_name)
    
    # Sort the models

    model_list = list(set(model_list))  # Remove duplicates if any
    sorted_models = sorted(model_list, key=extract_info)
    return sorted_models