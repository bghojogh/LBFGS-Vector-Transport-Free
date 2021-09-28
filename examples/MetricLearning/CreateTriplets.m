function [Triplets] = CreateTriplets(X, Y, tripletsize_per_class)

    UY = unique(Y);
 
    [ nx , dx ] = size(X);
    n_classes = length(UY);
    indices_of_points_in_classes = cell(n_classes,1);
    for class_ind=1:n_classes
        class_label = UY(class_ind);
        indices_of_points_in_classes{class_ind} = find(Y==class_label);
    end
    
    Triplets = zeros(tripletsize_per_class*n_classes, 3);
    triplet_index = 1;
    for class_ind=1:n_classes
        start_triplet_index = triplet_index;
        if tripletsize_per_class > length(indices_of_points_in_classes{class_ind})
            error("Error: The triplet size per class is larger than the class size!")
        end
        %%%% anchor-positive pairs:
        indices_in_this_class = indices_of_points_in_classes{class_ind};
        anchor_indices_in_dataset = randsample(indices_in_this_class, tripletsize_per_class);
        for anchor_index_in_dataset = anchor_indices_in_dataset'
            possible_positive_indices = indices_in_this_class(indices_in_this_class~=anchor_index_in_dataset);
            positive_index_in_dataset = randsample(possible_positive_indices, 1);  %--> sample without replacement
            Triplets(triplet_index, :) = [anchor_index_in_dataset, positive_index_in_dataset, 0];
            triplet_index = triplet_index + 1;
        end
        end_triplet_index = triplet_index - 1;
        %%%% negative samples:
        class_label = UY(class_ind);
        indices_of_points_in_enemy_classes = find(Y~=class_label);
        Triplets(start_triplet_index:end_triplet_index, 3) = randsample(indices_of_points_in_enemy_classes, tripletsize_per_class);
    end
    
end