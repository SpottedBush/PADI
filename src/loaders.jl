# 
# loaders.jl
# 
# Loading everything needed for the reconstruction. (TransTables etc...)
#----------------------------------------------------------
#
# This file is part of PADI
#
#
# Copyright (c) 2024-2025 Vincent Tardieux and Laurence Denneulin (see LICENCE.md)
#
#------------------------------------------------

const Trans_Table = Vector{NTuple{2, AffineTransform2D}}(); # Contains all affine transforms used 
const Star_Disk_Table = Vector{NTuple{4, AffineTransform2D}}(); # Contains all affine transforms used 
const Parameters = parameters_table[];
get_par()::parameters_table = Parameters[1];
const dataset = data_table[];


const EPSILON_save = Array{Float64, 1}(undef, 1);
function set_epsilon(epsilon)
    EPSILON_save[1]=epsilon
end
get_epsilon()::Float64=EPSILON_save[1];
const PRECOND_SAVE= Vector{Any}();
U()=PRECOND_SAVE[1];

const MASK_save = Vector{Array{Float64, 2}}();
get_MASK()=MASK_save[1];

function Indices(S::Int64,Nrot::Int64,Nframe::Int64)
#S = le nombre total de frame gauche ou droite
#Nrot : Nombre de positions du derotateur
#Nframe : Nombre de frame par positions de la lame demi-onde
	
	Nangle=Int32.(S/(Nrot*Nframe*4)) #Nombre de rotations de la lame demi-onde
	INDICES=zeros(4,Nframe*Nangle*Nrot)
	for i=1:4
		ind=repeat(range(0,stop=4*Nframe*(Nrot*Nangle-1),length=Nrot*Nangle), inner=Nframe)+(Nframe*i .-mod.(range(1,stop=Nframe*Nangle*Nrot,length=Nframe*Nangle*Nrot),Nframe))
		INDICES[i,:]=ind
	end
	INDICES=round.(Int64, INDICES);
	return INDICES
end

function trans_rotate(A::AffineTransform2D{Float64}, EPSILON, CENTER, NEWCENTER; ANGLE=0, is_star=true::Bool)
    if is_star
        return translate(translate(CENTER[1]-EPSILON[1],CENTER[2]-EPSILON[2], A),-NEWCENTER[1], -NEWCENTER[2])
    end
	return translate(rotate(translate(CENTER[1]-EPSILON[1],CENTER[2]-EPSILON[2], A), ANGLE),-NEWCENTER[1], -NEWCENTER[2])
end

function get_max_boxing(Id::AffineTransform2D{Float64}, size_data::NTuple{3, Int64}, center::Array{Float64, 1}, epsilon::Vector{NTuple{2, Array{Float64, 1}}}; derotang=nothing::Vector{Float64})
    bbox_output = (0, 0)
    if derotang !== nothing
        for k=1:size_data[3]
            A_left_star=trans_rotate(Id, epsilon[k][1], center, center)
            A_left_disk=trans_rotate(Id, epsilon[k][1], center, center, ANGLE=derotang[k], is_star=false)
            A_right_star=trans_rotate(Id, epsilon[k][2], center, center)
            A_right_disk=trans_rotate(Id, epsilon[k][2], center, center, ANGLE=derotang[k], is_star=false)

            out_left_star=bbox_size((size_data[1], Int64(size_data[2] / 2)), A_left_star)[1]
            out_left_disk=bbox_size((size_data[1], Int64(size_data[2] / 2)), A_left_disk)[1]
            out_right_star=bbox_size((size_data[1], Int64(size_data[2] / 2)), A_right_star)[1]
            out_right_disk=bbox_size((size_data[1], Int64(size_data[2] / 2)), A_right_disk)[1]

            xmax_star=max(bbox_output[1], out_left_star[1], out_right_star[1]);
            xmax_disk=max(bbox_output[1], out_left_disk[1], out_right_disk[1]);
            ymax_star=max(bbox_output[2], out_left_star[2], out_right_star[2]);
            ymax_disk=max(bbox_output[2], out_right_disk[2], out_right_disk[2]);

            bbox_output = (max(xmax_star, xmax_disk), max(ymax_star, ymax_disk));
        end
    else
        for k=1:size_data[3]
            A_left=translate(epsilon[k][1][1], epsilon[k][1][2], Id)
            A_right=translate(epsilon[k][2][1], epsilon[k][2][2], Id)
    
            out_left=bbox_size((size_data[1], Int64(size_data[2] / 2)), A_left)[1]
            out_right=bbox_size((size_data[1], Int64(size_data[2] / 2)), A_right)[1]
    
            xmax=max(bbox_output[1], out_left[1], out_right[1]);
            ymax=max(bbox_output[2], out_left[2], out_right[2]);
    
            bbox_output = (xmax,ymax);  
        end
    end
    return bbox_output

end

function push_to_trans_table(Id::AffineTransform2D{Float64}, epsilon::Vector{NTuple{2, Array{Float64, 1}}}, center::Array{Float64, 1}, newcenter::Tuple{Float64, Float64}; derotang=nothing::Vector{Float64})
    @assert isnothing(derotang) || length(epsilon) == length(derotang)
    n_frames = length(epsilon)
    if isnothing(derotang)
        for k=1:n_frames
            A_left=translate(epsilon[k][1][1] + centerdiff[1], epsilon[k][1][2] + centerdiff[2], Id)
            A_right=translate(epsilon[k][2][1] + centerdiff[1], epsilon[k][2][2] + centerdiff[2], Id)  
            push!(Trans_Table, (A_left, A_right))
        end
    else
        for k=1:n_frames
            A_left=inv(TransRotate(Id, epsilon[k][1], derotang[k], center, newcenter))
            A_right=inv(TransRotate(Id, epsilon[k][2], derotang[k], center, newcenter))   
            push!(Trans_Table, (A_left, A_right))
        end
    end
end

function push_to_star_disk_table(Id::AffineTransform2D{Float64}, epsilon::Vector{NTuple{2, Array{Float64, 1}}}, center::Array{Float64, 1}, newcenter::Array{Float64, 1}, derotang::Vector{Float64})
    for k=1:length(epsilon)
        A_left_star=inv(trans_rotate(Id, epsilon[k][1], center, newcenter))
        A_left_disk=inv(trans_rotate(Id, epsilon[k][1], center, newcenter, ANGLE=derotang[k], is_star=false))
        A_right_star=inv(trans_rotate(Id, epsilon[k][2], center, newcenter))
        A_right_disk=inv(trans_rotate(Id, epsilon[k][2], center, newcenter, ANGLE=derotang[k], is_star=false))
        push!(Star_Disk_Table, (A_left_star, A_left_disk, A_right_star, A_right_disk))
    end
end

function load_parameters(size_data::NTuple{3,Int64}, 
        Nframe::Int64,
        Nrot::Int64, 
        Nangle::Int64, 
        center::Array{Float64,1}, 
        psf_center::NTuple{2,Array{Float64,1}}, 
        epsilon::Vector{NTuple{2,Array{Float64,1}}};
        derotang=nothing::Vector{Float64},
        size_object=(0, 0)::NTuple{2,Int64},
        padding=20::Int64,
        v=nothing::Union{Array{Float64,2},Nothing})
    Id = AffineTransform2D{Float64}()
    sets_indices=Indices(size_data[3], Nangle, Nframe)
    if v !== nothing
        sets_v=Set_Vi(Nframe, size_data[3], v)
    else
        sets_v=Set_Vi(sets_indices)
    end
    bbox_output = get_max_boxing(Id, size_data, center, epsilon, derotang=derotang)
    if  (bbox_output[1] + padding > size_object[1]) || (bbox_output[2] + padding > size_object[2])
        @warn "The reconstruction size estimated by the mehod is larger than the size given by the user." "bbox_ouput set to the size obtained by the method." 
    end
    bbox_output = max(bbox_output .+ padding, size_object); 
    push!(Parameters, parameters_table((bbox_output[1], bbox_output[2], 4),
                    (size_data[1], size_data[2]), size_data[3], Nframe, Nrot, Nangle,
                    sets_v, sets_indices, center, psf_center, epsilon, derotang));
    newcenter = (bbox_output .+1)./2
    if size_data[3] == Nrot
        push_to_star_disk_table(Id, epsilon, center, [newcenter[1], newcenter[2]], derotang) 
        SetCropOperator(Star_Disk_Table)

    else
        push_to_trans_table(Id, epsilon, center, newcenter, derotang=derotang)
        SetCropOperator(Trans_Table)
    end
end

function load_data(name_data, name_weight)
    data=readfits(name_data);
    weight=readfits(name_weight);
    ker= MyKer;
    for k=1:size(data)[3]
        output_size = (get_par().rows[1], Int64(get_par().rows[2]/2));
        input_size = get_par().cols[1:2];
        T_l_star = TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, Star_Disk_Table[k][1])
        T_l_disk = TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, Star_Disk_Table[k][2])
        T_r_star = TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, Star_Disk_Table[k][3])
        T_r_disk = TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, Star_Disk_Table[k][4])
        push!(dataset, data_table((data[:,:,k]')[:,:], 
                    (weight[:,:,k]')[:,:],
                    FieldTransformOperator(get_par().cols,
                                            get_par().rows,
                                            get_par().v[k][1],
                                            get_par().v[k][2],
                                            T_l_star,
                                            T_l_disk,
                                            T_r_star,
                                            T_r_disk)));
    end
end

function flush_dataset()
    empty!(dataset)
end