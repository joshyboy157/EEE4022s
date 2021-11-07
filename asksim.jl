using Plots
plotly()
using Distributions
using FFTW

input = rand(0:1,10000)

#first block

f = 40000
sample_rate = 100000
A = 1
Δt = 1/sample_rate
chip_duration = 0.001 
nbits = length(input)
nsamples = Int(round(nbits*chip_duration/Δt))
y = zeros(nsamples)

#@show size(y)


function waveform(t)
    bit_num = Int(floor(t/chip_duration)+1)
    if(bit_num>nbits)
        return 0
    end
    b = input[bit_num]
    if b == 1
        y = cos(2*pi*f*t)
    else
        y = 0
    end
    return y
end


for n = 1:nsamples
    t = (n-1)*Δt
    y[n] = waveform(t)
end

#next block
#display(plot(y[1:2000]))

# #next block

Y = fft(y)
n = length(Y)
f_axis = (0:n-1)

# #next block

#display(plot(f_axis,abs.(Y)))

# # next block bandlimit
H = zeros(n)
H[385000:415000] .=1
H[585000:615000] .=1

#display(plot(f_axis,H))

#next block

noise = zeros(nsamples)
Sigma = [0.7, 1, 2, 3, 4, 5, 7,9 ,13, 15, 20 ,30 ,40,50,70,100,130,150,180,200,500,1000,2000]
noise1=[zeros(nsamples),zeros(nsamples),zeros(nsamples)]
sigma = 150*A

for u in 1:length(Sigma)
    noise = zeros(nsamples)
    for n = 1:nsamples
        noise[n] = Sigma[u]*rand(Uniform(-1,1))
    end
    pushfirst!(noise1,noise)
end

pop!(noise1)
pop!(noise1)
pop!(noise1)
# display(plot(noise1[1]))
# display(plot(noise1[10]))
# display(plot(noise1[19]))
#display(plot(noise))
noise1 = reverse(noise1)
b= fft(noise1[1])

# #next block
Noise1 = [b]
for i in 2:length(noise1)
    pushfirst!(Noise1,fft(noise1[i]))
end


Noise1 = reverse(Noise1)


# #next block


# #next block


lim_Noise1 = [H.*abs.(Noise1[1])]
for i in 2:length(noise1)
    pushfirst!(lim_Noise1,H.*abs.(Noise1[i]))
end

lim_Noise1 = reverse(lim_Noise1)



# #next block

# display(plot(f_axis,abs.(lim_Noise1[1]),xlabel = "Frequency Hz" , ylabel = "Magnitude"))
# display(plot(f_axis,abs.(lim_Noise1[13])))
# display(plot(f_axis,abs.(lim_Noise1[19])))


# #next block

lim_noise1 = [abs.(ifft(lim_Noise1[1]))]
#pushfirst!(lim_noise1,abs.(ifft(lim_Noise1[2])))

for i in 2:length(noise1)
    b = abs.(ifft(lim_Noise1[i]))
    pushfirst!(lim_noise1,b)
end

lim_noise1 = reverse(lim_noise1)

# display(plot(lim_noise1[1][4000:200000]))
# display(plot(lim_noise1[10][4000:200000]))
# display(plot(lim_noise1[15][4000:200000]))
# display(plot(lim_noise1[20][4000:200000]))
# #next block

sent = [y.+ lim_noise1[1]]
for i in 2:length(noise1)
    pushfirst!(sent,y.+ lim_noise1[i])
end

sent = reverse(sent)

# display(plot(sent[2][4000:5000]))
# display(plot(sent[9][4000:5000]))
# display(plot(sent[10][4000:5000]))
# display(plot(sent[20][4000:5000]))
# #next block

SENT = fft(sent[1])

# #next block

#display(plot(sent[3000:900000]))


# #next block

sizeofchip = nsamples/nbits
sizeofchip = Int(round(sizeofchip))
its = nbits
demod = zeros(nbits)
demod = complex(demod)
BER_arr = [1.0]
pop!(BER_arr)
std_y = [1.0]
std_noise = [1.0]
pop!(std_y)
pop!(std_noise)

for i in 1:length(sent)

    for num in 1:its
        part = sent[i][(num-1)*sizeofchip+1:num*sizeofchip]
        t = ((num-1)*sizeofchip:num*sizeofchip-1)*Δt
        iq = sum(part .* exp.(2*π*f*t*im))*Δt
        demod[num] = iq
    end

    # @show length(demod)
    # @show demod[5]
    # @show abs(demod[5])
    # @show demod[51]
    # @show abs(demod[51])
    # @show demod[511]
    # @show abs(demod[511])
    # #next block

    output_bits = zeros(nbits)
    for i in 1:length(demod)
    #     println(abs(i))
     if  abs(demod[i])>0.00025
            output_bits[i] = 1
        else
            output_bits[i] = 0
        end
    end

    #@show output_bits[5]
    #@show input[5]

    # #next block

    BER = 0
    for n = 1:nbits
        if input[n]!= output_bits[n]
            BER+=1
        end
    end

    BER = BER/nbits
    pushfirst!(BER_arr,BER)
    pushfirst!(std_y,std(y))
    pushfirst!(std_noise,std(lim_noise1[i]))
    # # #next block
    # #print("power of y: ")
    # println(std(y))
    # #print("power of noise: ")
    # println(std(abs.(lim_noise)))
    # println(BER)



    # for i in 1:27
    #     b = std_y[i]
    #     e = std_noise[i]
    #     t = b/e
    #     t=t*t
    #     pushfirst!(SNR,t)
    # end

    # pop!(SNR)
    # println(length(SNR))
    # SNR = reverse(SNR)
    # SNRDec = log10.(SNR)
    # print(SNRDec)


    # display(plot(SNRDec,BER_arr,xlabel ="SNR(Db)" , ylabel = "BER"))
end

std_y = reverse(std_y)
std_noise = reverse(std_noise)
BER_arr = reverse(BER_arr)
println(BER_arr)
println(std_y)
println(std_noise)
SNR = [std_y[1]/std_noise[1]]

for i in 2:length(std_y)
    o = std_y[i]/std_noise[i]
    pushfirst!(SNR,o)
end

SNR = reverse(SNR)

SNRdB = log10.(SNR)

display(plot(SNRdB,BER_arr,xlabel = "SNR dB" , ylabel = "BER"))

