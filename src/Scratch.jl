using CSV

infile = "D:\\Dropbox\\dev\\SFTSV\\SFTSV_S_info.csv"
fin = open(infile, "r")
lines = readlines(fin)
header = split(lines[1],",")
age = Int[]
death = Int[]
names = String[]
sequences = String[]
for i=2:length(lines)
    if length(lines[i]) > 0
        spl = split(lines[i],",")
        if length(spl) >= 14
            if length(spl[3]) > 0
                push!(age, parse(Int, spl[3]))
                push!(death, parse(Int, spl[8]))
                push!(names, spl[11])
                push!(sequences, spl[14])
            end
        end
    end
end

allages = true
outfile = string(infile,".under60.fas")
if allages
    outfile = string(infile,".allages.fas")
end

fout = open(outfile, "w")
csvout = open(string(outfile,".tsv"),"w")
for i=1:length(names)
    deathlabel = "S"
    if death[i] == 1
        deathlabel = "D"
    end
    #=
    deathage = "ND60"
    if death[i] == 1 && age[i] < 60
        deathage = "DU60"
    end=#
    if allages || age[i] < 60
        newname = string(names[i][2:end],"_",deathlabel,"_",age[i])
        write(fout, string(">",newname,"\n"))
        write(fout, sequences[i],"\n")
        if death[i] == 1
            write(csvout, string(newname,"\t1\n"))
        else
            write(csvout, string(newname,"\t0\n"))
        end
    end
end
close(csvout)
close(fout)
