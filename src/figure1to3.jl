using CairoMakie
using FileIO
include("error_ellipse.jl")
include("therapy_ebs.jl")

fs = 30 # fontsize
fn = "Calibri"

r = rand(d_obs,10000)

# ==================================================
f1 = Figure(resolution=(855,422),fontsize = fs,font=fn)
# Column 1, Row 1:
f1[1,1] = Axis(f1,title = "log(HDL-Ch)", limits=(-1, 1.5, 0, 2),
    yticks=[0,1,2],xticks=[-1,0,1],ylabel="PDF")

sim = hist!(f1[1,1],obs_pp[1,selidx], normalization = :pdf) # Simulation
#hist!(f1[1,1][1,1],randn(1000), normalization = :pdf) # Data
da = density!(f1[1,1],r[1,:],color=(:red,0),strokecolor=:red,strokewidth=3)

# Column 1, Row 2:
f1[2,1] = Axis(f1,xlabel="log(HDL-Ch)",ylabel="log(TPCh)", limits = (-1, 1.5, 0.5, 2.5),yticks=[0,1,2,3],xticks=[-1,0,1])
scatter!(f1[2,1], obs_pp[1,selidx],obs_pp[2,selidx]) # Scatter of simulation
d = return_nhanes([1;2])
(x,y) = getellipsepoints(d.μ,d.Σ,0.95)
lines!(f1[2,1],x,y,linewidth=4,color=:red)

# Column 2, Row 1:
f1[1,2] = Axis(f1, title="log(TPCh)", limits=(0.5, 2.5, 0, 2),yticks=[0,1,2],xticks=[0,1,2,3])
hist!(f1[1,2],obs_pp[2,selidx], normalization = :pdf) # Simulation
density!(f1[1,2],r[2,:],color=(:red,0),strokecolor=:red,strokewidth=3)
#hist!(f1[1,2][1,1],randn(1000), normalization = :pdf) # Data
# Column 2, Row 2:
f1[2,2] = Axis(f1, xlabel = "log(TPCh)", ylabel="log(TPTG)", limits = (0.5, 2.5, -2, 2),xticks=[0,1,2,3],yticks=[-2,0,2])
scatter!(f1[2,2], obs_pp[2,selidx],obs_pp[3,selidx]) # Scatter of simulation
d = return_nhanes([2;3])
(x,y) = getellipsepoints(d.μ,d.Σ,0.95)
lines!(f1[2,2],x,y,linewidth=4,color=:red)

# Column 3, Row 1:
f1[1,3] = Axis(f1,title="log(TPTG)", limits=(-2, 2, 0, 2),yticks=[0,1,2],xticks=[-2,0,2])
hist!(f1[1,3],obs_pp[3,selidx], normalization = :pdf) # Simulation
density!(f1[1,3],r[3,:],color=(:red,0),strokecolor=:red,strokewidth=3)
#hist!(f1[1,3][1,1],randn(1000), normalization = :pdf) # Data
# Column 3, Row 2:
f1[2,3] = Axis(f1,ylabel="log(HDL-Ch)", xlabel="log(TPTG)", limits = (-2, 2, -1, 1.5),xticks=[-2,0,2],yticks=[-1,0,1])
scatter!(f1[2,3], obs_pp[3,selidx],obs_pp[1,selidx]) # Scatter of simulation
d = return_nhanes([1;3])
(x,y) = getellipsepoints(d.μ,d.Σ,0.95)
lines!(f1[2,3],y,x,linewidth=4,color=:red) # Note the reverse y, x

#f1[3,1] = Axis(f1)
elem_1 = [LineElement(color = :red, linestyle = nothing,linewidth=5)]
Legend(f1[3,2],[sim,elem_1],["Virtual Patients","NHANES 2018 [6]"],orientation=:horizontal,
    tellheight=true, tellwidth = false, framevisible = false)

# Format, display, and save:
#resize_to_layout!(f1)
f1
FileIO.save("./fig/figure1.png",f1)
# ==================================================

# ==================================================
f2 = Figure(resolution=(855,422),fontsize = fs,font=fn)
# Col 1, liver fat scatter:
f2[1,1] = Axis(f2,title= "Liver TG", xlabel="LTG [% of Liver]", ylabel="PDF",
    yticks = [0,0.1,0.2],limits = (0,25,0,0.2))
hist!(f2[1,1],exp.(obs_pp[4,selidx]),normalization=:pdf)
density!(f2[1,1],exp.(r[4,:]),color=(:red,0),strokecolor=:red,strokewidth=3)
poly!(f2[1,1],Point2f[(5,0),(25,0),(25,0.2),(5,0.2)],color=(:red,0.1))
text!(f2[1,1], 9.5, 0.18,text="NAFLD",textsize=30,align=(:center,:center))

# Col 2, liver fat vs. TG scatter:
f2[1,2] = Axis(f2,title= "Liver vs. Plasma TG", ylabel="log(TPTG)", xlabel="log(LTG)", limits=(-2,5,-2,2),xticks=[-2,0,2,4],yticks=[-2,0,2])
scatter!(f2[1,2],(obs_pp[4,selidx]),(obs_pp[3,selidx]))
d = return_nhanes([3;4])
(x,y) = getellipsepoints(d.μ,d.Σ,0.95)
lines!(f2[1,2],(y),(x),linewidth=4,color=:red) # Note the reverse y, x

Legend(f2[2,1:2],[sim,elem_1],["Virtual Patients","Data [6,7]"],orientation=:horizontal,
    tellheight=true, tellwidth = false, framevisible = false)


# Format, display, and save:
resize_to_layout!(f2)
f2
FileIO.save("./fig/figure2.png",f2)
# ==================================================

# ==================================================
colors = Makie.wong_colors()

ebs = return_therapy_ebs()

x = [1,2,3,1,2,3] # X-axis position
y_cetp = [mean(hdl_pcb_cetpi[selidx]);mean(tc_pcb_cetpi[selidx]);mean(tg_pcb_cetpi[selidx]);22;0.48;-4] # Simulated/Data
y_pcsk9 = [mean(hdl_pcb_pk9ab[selidx]);mean(tc_pcb_pk9ab[selidx]);mean(tg_pcb_pk9ab[selidx]);7;-40;-17.4]
grp = [1,1,1,2,2,2] # Distringuish simulated from data

xVP = [1-0.21;1.8;3.0-0.22]
xData = [1.22;2.2;3.22]

f3 = Figure(resolution=(855,422),fontsize = fs,
            font=fn)
# Col 1, dalcetrapib:
f3[1,1] = Axis(f3,title="Dalcetrapib\n600mg QD", 
            ylabel="Δ%Pbo, 4wks",
            xticks = (1:3,["HDL-Ch","TPCh","TPTG"]),limits=(0.5,3.5,-60,40))
barplot!(f3[1,1],xData,y_cetp[4:6],color=:red, width = 0.45)
#barplot!(f3[1,1],x,y_cetp,
#    dodge = grp,
#    colormap = [colors[1],:red],
#    color = grp)

nvp = sum(selidx)
boxplot!(f3[1,1],repeat([xVP[1]],nvp),hdl_pcb_cetpi[selidx],width = 0.475,color=colors[1])
boxplot!(f3[1,1],repeat([xVP[2]],nvp),tc_pcb_cetpi[selidx],width = 0.475,color=colors[1])
boxplot!(f3[1,1],repeat([xVP[3]],nvp),tg_pcb_cetpi[selidx],width = 0.475,color=colors[1])

# Add Data EBs:
errorbars!(f3[1,1],xData,y_cetp[4:6],y_cetp[4:6] .- ebs[1,1:3],ebs[2,1:3] .- y_cetp[4:6], whiskerwidth = 10, linewidth = 2)

#barplot!(f3[1,1],[1.5,3.5,5.5], randn(3), width = 1)
# Col 2, Evolocumab:
f3[1,2] = Axis(f3,title="Evolocumab\n140mg Q2W", 
            ylabel="Δ%Pbo, 12wks",
            xticks = (1:3,["HDL-Ch","TPCh","TPTG"]),limits=(0.5,3.5,-60,40))

barplot!(f3[1,2],xData,y_pcsk9[4:6],color=:red, width = 0.45)
#=
barplot!(f3[1,2], x, y_pcsk9,
    dodge = grp,
    colormap = [colors[1],:red],
    color = grp)
    =#
errorbars!(f3[1,2],xData,y_pcsk9[4:6],y_pcsk9[4:6] .- ebs[1,4:6],ebs[2,4:6] .- y_pcsk9[4:6], whiskerwidth = 10, linewidth = 2)

boxplot!(f3[1,2], repeat([xVP[1]],nvp), hdl_pcb_pk9ab[selidx],width = 0.475,color=colors[1])
boxplot!(f3[1,2],repeat([xVP[2]],nvp),tc_pcb_pk9ab[selidx],width = 0.475,color=colors[1])
boxplot!(f3[1,2],repeat([xVP[3]],nvp),tg_pcb_pk9ab[selidx],width = 0.475,color=colors[1])


elem_1 = [PolyElement(polycolor = colors[1])]
elem_2 = [PolyElement(polycolor = :red)]
Legend(f3[2,1:2],[elem_1, elem_2],["VPs","Data, Mean±95%CIs [8,9]"],orientation=:horizontal,
    tellheight=true, tellwidth = false, framevisible = false)
# Format, display, and save:
resize_to_layout!(f3)
f3
FileIO.save("./fig/figure3.png",f3)
# ==================================================
f4 = Figure(resolution=(855,422),fontsize = fs,font=fn)
f4[1,1] = Axis(f4,title="Sensitivity At Baseline",
    xticks=([1,2,3,4],["HDL-Ch","TPCh","TPTG","LTG"]), 
    yticks = (1:sum(scal_r),sdSA[scal_r]), xlabel="End Point")
hm = heatmap!(f4[1,1], hm_scal[scal_r,:]',colorrange=(0,1),colormap=:temperaturemap)
Colorbar(f4[1,2],hm,ticks=([0,0.5,1.0],["(-)", "(0)", "(+)"]))

resize_to_layout!(f4)
f4
FileIO.save("./fig/figure4.png",f4)