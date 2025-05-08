using PADI
using DelimitedFiles
using EasyFITS
using InterpolationKernels
using Statistics: mean

function Double_Difference(data::Array{Float64,4},ind::Array{Int64,2})
    n1,n2,n3,n4= size(data);
    S=Array{PolarimetricPixel{Float64}, 2}(undef,n1,n2);
    @inbounds for i2 in 1:n2
        @simd for i1 in 1:n1   
            S[i1,i2]=Double_Difference(data[i1,i2,:,:], ind);
        end
    end
    return PolarimetricMap(S)
end

function Double_Difference(data::Array{Float64,2},ind::Array{Int64,2})
    Q=mean((data[ind[1,:],1] .-data[ind[1,:],2] 
              .-(data[ind[2,:],1] .-data[ind[2,:],2]))/2);		
    U=mean((data[ind[3,:],1] .-data[ind[3,:],2] 
              .-(data[ind[4,:],1] .-data[ind[4,:],2]))/2);
    Iq=(data[ind[1,:],1] .+data[ind[1,:],2] 
        .+data[ind[2,:],1] .+data[ind[2,:],2])/2;
    Iu=(data[ind[3,:],1] .+data[ind[3,:],2] 
        .+data[ind[4,:],1] .+data[ind[4,:],2])/2;

    return PolarimetricPixel("stokes", 
                             (mean(Iq)+ mean(Iu))/2,
                             0.0, 
                             Q,
                             U)
end


contrast_list = [-2.0]
# Make sure the folder exists
for k in contrast_list
    mkpath("test_results/contrast_10e$(k)/")
end

par=readdlm("data_for_demo/Parameters.txt")
DSIZE=Int64(par[1]);
NTOT=Int64(par[2]);
Nframe=Int64(par[3]);
Nrot=Int64(par[4]);
Nangle=NTOT÷(Nframe*4)
Center=par[5:6];
# DerotAng = deg2rad.(readdlm("data_for_demo/pds70_angles.txt", Float64)[1:64])
max_angle = 64
DerotAng = [deg2rad(i) for i in range(1, max_angle, length=64)]
Epsilon=Vector{Tuple{Array{Float64,1},Array{Float64,1}}}();

for iter=1:NTOT
    ind=div(iter-1, NTOT/4)+1
    push!(Epsilon,([0. ,0. ],par[end-1:end]));
end

psf_center=readdlm("data_for_demo/PSF_centers_Airy.txt");

PADI.load_parameters((DSIZE, 2*DSIZE, NTOT), Nframe, Nrot
, Nangle, Center, (psf_center[1:2], psf_center[3:4]), Epsilon, derotang=DerotAng)

for k in contrast_list
    load_data("test_results/contrast_10e$(k)/DATA.fits", 
            "test_results/contrast_10e$(k)/WEIGHT.fits")
	Sdim=length(PADI.dataset)
    
    DATA=zeros(get_par().cols[1], get_par().cols[2], Sdim,2);
    WEIGHT=zeros(get_par().cols[1], get_par().cols[2], Sdim,2);

    ker= CatmullRomSpline(Float64, Flat)       
    input_size=(get_par().rows[1], get_par().rows[2]÷2);
    output_size= get_par().cols[1:2];

    # Pre processing
    for i=1:Sdim
        T1=TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, inv(PADI.Star_Disk_Table[i][2]))
        T2=TwoDimensionalTransformInterpolator(output_size, input_size, ker, ker, inv(PADI.Star_Disk_Table[i][4]))

        data = PADI.dataset[i].data
        weights = PADI.dataset[i].weights 

            for i=2:size(data)[1]-1
                for j=2:size(data)[2]-1
                    if weights[i,j] == 0.
                        data[i,j] = (data[i-1, j] +data[i+1, j] +data[i, j-1] + data[i, j+1]) /
                                ((data[i-1, j] !=0) + (data[i+1, j]!=0) + (data[i, j-1]!=0) + (data[i, j+1]!=0))
                    end
                end
             end

        I1=T1 * data[:,1:end÷2]
        I2=T2 * data[:,1+end÷2:end]

        DATA[:,:,i,1]=I1;
        DATA[:,:,i,2]=I2;
    end
    DD=Double_Difference(DATA, PADI.Indices(NTOT, 1, Nframe));
    write_polar_map(DD, "test_results/contrast_10e-2.0/Results_Separable_DoubleDifference.fits", overwrite=true)
end