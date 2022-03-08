#---------------------Packages needed to run this script------------------
# install.packages("ncdf4")
#----------------------------------------------------------
# This script requires a netcdf trajectory of the oxygen atoms of water molecules and 12 reference atoms that will be used to trace 4 cylinders in order to track the water molecules (i.e. oxygen atoms) in each pore in each frame. This script was written to process permeabilities of each monomer in a tetrameric arrangement, so every calculation is done 4 times. Permeabilities are calculated according "Collective diffusion model for water permeation through microscopic channels" by Zhu, Tajkhorshid and Schulten, 2004 #

library(ncdf4) # package for netcdf manipulation

# while (TRUE) {
#   ncfname <- readline(prompt = "Enter file name: ")
#   if (!(ncfname %in% list.files())) {
#     print("Invalid filename")
#   } else {
#     break()
#   }
# }
#
# while (TRUE) {
#   monomer_count <-
#     as.integer(readline(prompt = "How many monomers are in the structure?: "))
#   if (is.na(monomer_count)) {
#     print("Invalid input. Enter an integer number")
#   } else {
#     print(paste("The amount of monomers are", monomer_count))
#     break()
#   }
# }


ncfname <- paste("./AOX.nc")
monomer_count <- 4
ref_atoms_per_monomer <- 3
ref_atoms_names <- c(1:ref_atoms_per_monomer)
monomer_names <- toupper(letters[1:monomer_count])
ncin <- nc_open(ncfname)
data_name <- "coordinates"
atoms_length <- ncin[["var"]][["coordinates"]][["size"]][2]
time_length <- ncin[["var"]][["coordinates"]][["size"]][3]

# Reference atoms naming
# reference atoms must be the first in the list, just before the oxygen atoms. If not, change the index number in the variables below
reference_1 <- 1
reference_2 <- 2
reference_3 <- 3
frame_step <- 1 # in picoseconds
iii <- -1


for (monomer in monomer_names) {
  iii <- iii + 1
  for (ref_atom in ref_atoms_names) {
    assign(paste("ref_", ref_atom, monomer, sep = ""), get(paste("reference_", ref_atom, sep = "")) + 3 * iii)
  }
  assign(
    paste("WAT_", monomer, "_atoms", sep = ""),
    data.frame(matrix(ncol = 0, nrow = 0))
  )
}

# this will find any oxygen atom that appears inside the cylinder defined by the reference atoms and a constant radius. Note that each cylinder center and length is defined for each frame, so it has tiny length fluctuations over time, but with this we can be sure that the cylinder follows the pore along the trajectory.

pedazos <-
  20 # number of fragments (in spanish: "pedazos") in which the trajectory will be cut for this analysis. This number depends on the RAM of your system, the length of the trajectory and the number of oxygen atoms in your system. You will have to benchmark this. Start with 50 and check your RAM usage. If it is low, then you can explore lower values. Lower values generates faster calculations and bigger RAM usages.

print("buscando_agua_dentro_del_poro")
for (nnn in (0:(pedazos - 1))) {
  print(paste("procesando_pedacito", nnn + 1, "de", pedazos))

  xyz <-
    ncvar_get(
      ncin,
      "coordinates",
      start = c(1, 1, 1 + nnn * (time_length / pedazos)),
      count = c(3, -1, time_length / pedazos)
    )

  for (monomer in monomer_names) {
    assign(
      paste("z_1", monomer, sep = ""),
      ncvar_get(
        ncin,
        "coordinates",
        start = c(3, get(paste(
          "ref_1", monomer,
          sep = ""
        )), 1 + nnn * (time_length / pedazos)),
        count = c(1, 1, time_length / pedazos)
      )
    )

    assign(
      paste("xyz_2", monomer, sep = ""),
      ncvar_get(
        ncin,
        "coordinates",
        start = c(1, get(paste(
          "ref_2", monomer,
          sep = ""
        )), 1 + nnn * (time_length / pedazos)),
        count = c(3, 1, time_length / pedazos)
      )
    )

    assign(
      paste("xy_3", monomer, sep = ""),
      ncvar_get(
        ncin,
        "coordinates",
        start = c(1, get(paste(
          "ref_3", monomer,
          sep = ""
        )), 1 + nnn * (time_length / pedazos)),
        count = c(2, 1, time_length / pedazos)
      )
    )

    assign(paste("z_top", monomer, sep = ""), get(paste("z_1", monomer, sep = "")))
    assign(paste("x_2", monomer, sep = ""), get(paste("xyz_2", monomer, sep = ""))[1, ])
    assign(paste("y_2", monomer, sep = ""), get(paste("xyz_2", monomer, sep = ""))[2, ])
    assign(paste("z_bot", monomer, sep = ""), get(paste("xyz_2", monomer, sep = ""))[3, ])
    assign(paste("x_3", monomer, sep = ""), get(paste("xy_3", monomer, sep = ""))[1, ])
    assign(paste("y_3", monomer, sep = ""), get(paste("xy_3", monomer, sep = ""))[2, ])
    assign(paste("x_cent"), (get(paste(
      "x_2", monomer,
      sep = ""
    )) + get(paste(
      "x_3", monomer,
      sep = ""
    ))) / 2)
    assign(paste("y_cent"), (get(paste(
      "y_2", monomer,
      sep = ""
    )) + get(paste(
      "y_3", monomer,
      sep = ""
    ))) / 2)
    assign(paste("z_max"), max(get(paste(
      "z_top", monomer,
      sep = ""
    ))))
    assign(paste("z_min"), max(get(paste(
      "z_bot", monomer,
      sep = ""
    ))))

    assign(
      paste("WAT_", monomer, "_atoms", sep = ""),
      rbind(
        get(paste("WAT_", monomer, "_atoms", sep = "")),
        as.matrix(unique(arrayInd(
          which(t((t(xyz[1, , ]) - c(x_cent))^2) + t((t(xyz[2, , ]) - c(y_cent))^2) < 36.0 & xyz[3, , ] < z_max & xyz[3, , ] > z_min), dim(xyz[1, , ])
        )[, 1]))
      )
    )
    # WAT_A_atoms <- rbind(WAT_A_atoms,as.matrix(unique(arrayInd(which(t((t(xyz[1,,])-c(x_centA))^2)+t((t(xyz[2,,])-c(y_centA))^2)<36 & xyz[3,,] < z_maxA & xyz[3,,] > z_minA),dim(xyz[1,,]))[,1])))
  }
}



rm(xyz)
gc()
for (monomer in monomer_names) {
  assign(
    paste("WAT_", monomer, "_atoms", sep = ""),
    as.data.frame(unique(get(
      paste("WAT_", monomer, "_atoms", sep = "")
    )[get(paste("WAT_", monomer, "_atoms", sep = "")) > monomer_count * 3]))
  )

  assign(paste("df", monomer, sep = ""), data.frame(matrix(ncol = 5, nrow = 0)))
  tmp <- get(paste("df", monomer, sep = ""))
  colnames(tmp) <- c("time", "atom", "x", "y", "z")
  assign(paste("df", monomer, sep = ""), tmp)

  tmp <- get(paste("WAT_", monomer, "_atoms", sep = ""))
  colnames(tmp) <- "atoms"
  assign(paste("WAT_", monomer, "_atoms", sep = ""), tmp)

  tmp <- get(paste("df", monomer, sep = ""))
  colnames(tmp) <- c("time", "atom", "x", "y", "z")
  assign(paste("df", monomer, sep = ""), tmp)
}


# This will load the trajectory of the every atom found in the prior step in a dataframe

pedacitos <- 20
# number of fragments (in spanish: "pedacitos") in which the trajectory will be cut for this analysis. This number depends on the RAM of your system, the length of the trajectory and the number of oxygen atoms in your system. You will have to benchmark this. Start with 50 and check your RAM usage. If it is low, then you can explore lower values. Lower values generates faster calculations and bigger RAM usages.

print("extrayendo_atomos_de_la_trayectoria")

for (m in (0:(pedacitos - 1)))
{
  print(paste("procesando_pedacito", m + 1, "de", pedacitos))
  atom_xyz <-
    ncvar_get(
      ncin,
      "coordinates",
      start = c(1, 1, (1 + m * (
        time_length / pedacitos
      ))),
      count = c(3, atoms_length, time_length / pedacitos)
    )
  for (monomer in monomer_names) {
    for (atom in get(paste("WAT_", monomer, "_atoms", sep = ""))[, 1]) {
      assign(
        paste("data", monomer, sep = ""),
        cbind(c((
          1 + m * (time_length / pedacitos)
        ):((1 + m) * (time_length / pedacitos)
        )), c(rep(
          atom, (time_length / pedacitos)
        )), t(atom_xyz[, as.integer(atom), ]))
      )
      assign(
        paste("df", monomer, sep = ""),
        rbind(get(paste(
          "df", monomer,
          sep = ""
        )), get(paste(
          "data", monomer,
          sep = ""
        )))
      )
    }
  }
}


rm(atom_xyz); rm(list=ls(pattern="^data")); rm(list=ls(pattern="^WAT"))
rm(list=ls(pattern="^xy"))




# In this step dz (refer to Zhu et al. 2004) will be computed for every frame for every atom whenever they are inside the pore. Z1-Z0 will be computed as dz for frame 1 for an atom if the atom appears in the pore in any of the 2 frames. If the atom is outside the pore in the two consecutive frames, then dz=0.

for (chain in monomer_names)
{
  tmp <- get(paste("df", chain, sep = ""))
  colnames(tmp) <- c("time", "atom", "x", "y", "z")
  assign(paste("df", chain, sep = ""), tmp)
  df <- as.data.frame(get(paste("df", chain, sep = "")))
  dfdz <- data.frame(matrix(ncol = 3, nrow = 0))
  contador <- 0
  
  if (length(df[, 2]) > 0){
    
    for (k in unique(df[, 2]))
    {
      contador <- contador + 1
      print(paste(
        "procesando_atomo",
        k,
        chain,
        contador,
        "de",
        length(unique(df[, 2]))
      ))
      nuevo <- df[which(df$atom == k), ]
      diferencias_z <- as.numeric(as.character(df[which(df$atom == k), ][2:time_length, 5])) - as.numeric(as.character(df[which(df$atom == k), ][1:(time_length - 1), 5]))
      diferencias_z[time_length] <- 0
      diferencias_z <- as.data.frame(diferencias_z)
      for (j in 2:time_length)
      {
        z_top <-
          ncvar_get(ncin,
                    "coordinates",
                    start = c(3, get(paste(
                      "ref_1", chain,
                      sep = ""
                    )), j),
                    count = c(1, 1, 1)
          )
        x_2 <-
          ncvar_get(ncin,
                    "coordinates",
                    start = c(1, get(paste(
                      "ref_2", chain,
                      sep = ""
                    )), j),
                    count = c(1, 1, 1)
          )
        y_2 <-
          ncvar_get(ncin,
                    "coordinates",
                    start = c(2, get(paste(
                      "ref_2", chain,
                      sep = ""
                    )), j),
                    count = c(1, 1, 1)
          )
        z_bot <-
          ncvar_get(ncin,
                    "coordinates",
                    start = c(3, get(paste(
                      "ref_2", chain,
                      sep = ""
                    )), j),
                    count = c(1, 1, 1)
          )
        x_3 <-
          ncvar_get(ncin,
                    "coordinates",
                    start = c(1, get(paste(
                      "ref_3", chain,
                      sep = ""
                    )), j),
                    count = c(1, 1, 1)
          )
        y_3 <-
          ncvar_get(ncin,
                    "coordinates",
                    start = c(2, get(paste(
                      "ref_3", chain,
                      sep = ""
                    )), j),
                    count = c(1, 1, 1)
          )
        x_cent <- (x_2 + x_3) / 2
        y_cent <- (y_2 + y_3) / 2
        x_atom <- nuevo[j, 3]
        y_atom <- nuevo[j, 4]
        z_atom <- nuevo[j, 5]
        z_top_p <-
          ncvar_get(ncin,
                    "coordinates",
                    start = c(3, get(paste(
                      "ref_1", chain,
                      sep = ""
                    )), j - 1),
                    count = c(1, 1, 1)
          )
        x_2_p <-
          ncvar_get(ncin,
                    "coordinates",
                    start = c(1, get(paste(
                      "ref_2", chain,
                      sep = ""
                    )), j - 1),
                    count = c(1, 1, 1)
          )
        y_2_p <-
          ncvar_get(ncin,
                    "coordinates",
                    start = c(2, get(paste(
                      "ref_2", chain,
                      sep = ""
                    )), j - 1),
                    count = c(1, 1, 1)
          )
        z_bot_p <-
          ncvar_get(ncin,
                    "coordinates",
                    start = c(3, get(paste(
                      "ref_2", chain,
                      sep = ""
                    )), j - 1),
                    count = c(1, 1, 1)
          )
        x_3_p <-
          ncvar_get(ncin,
                    "coordinates",
                    start = c(1, get(paste(
                      "ref_3", chain,
                      sep = ""
                    )), j - 1),
                    count = c(1, 1, 1)
          )
        y_3_p <-
          ncvar_get(ncin,
                    "coordinates",
                    start = c(2, get(paste(
                      "ref_3", chain,
                      sep = ""
                    )), j - 1),
                    count = c(1, 1, 1)
          )
        x_cent_p <- (x_2_p + x_3_p) / 2
        y_cent_p <- (y_2_p + y_3_p) / 2
        x_atom_p <- nuevo[j - 1, 3]
        y_atom_p <- nuevo[j - 1, 4]
        z_atom_p <- nuevo[j - 1, 5]
        if ((((x_atom - x_cent)^2 + (y_atom - y_cent)^2 >= 36) ||
             (z_bot >= z_atom) ||
             (z_atom >= z_top)) &&
            (((x_atom_p - x_cent_p)^2 + (y_atom_p - y_cent_p)^2 >= 36) ||
             (z_bot_p >= z_atom_p) || (z_atom_p >= z_top_p))) {
          diferencias_z[j, 1] <- 0
        }
      }
      data <-
        cbind(c(1:time_length), c(rep(k, time_length)), diferencias_z[, 1])
      dfdz <- rbind(dfdz, data)
        assign(paste("dfdz", chain, sep = ""), dfdz)
    }
    
  } else {
    
    assign(paste("dfdz", chain, sep = ""), data.frame(0,0,0))
 
  }
  
  
  
  
  rm(dfdz)

  tmp <- get(paste("dfdz", chain, sep = ""))
  colnames(tmp) <- c("time", "atom", "dz")
  assign(paste("dfdz", chain, sep = ""), tmp)
}



for (monomer in monomer_names) {
  # dn (Zhu et al 2004) are computed for every frame and every monomer.
  assign(paste("dn", monomer, sep = ""), data.frame(matrix(ncol = 2, nrow = 0)))
  assign(paste("dif"), data.frame(matrix(ncol = 1, nrow = 0)))
  for (iii in 1:time_length) {
    assign(paste("dif"), rbind(dif, sum(get(
      paste("dfdz", monomer, sep = "")
    )[which(get(paste("dfdz", monomer, sep = ""))$time == iii), 3])))
  }
  assign(paste("dn", monomer, sep = ""), cbind(c(1:time_length), dif))

  tmp <- get(paste("dn", monomer, sep = ""))
  colnames(tmp) <- c("time", "dn")
  assign(paste("dn", monomer, sep = ""), tmp)

  # L (pore length) is computed for every frame and each monomer
  assign(
    paste("longitud_tubo_", monomer, sep = ""),
    mean(
      ncvar_get(
        ncin,
        "coordinates",
        start = c(3, get(paste(
          "ref_1", monomer,
          sep = ""
        )), 1),
        count = c(1, 1, time_length)
      ) - ncvar_get(
        ncin,
        "coordinates",
        start = c(3, get(paste(
          "ref_2", monomer,
          sep = ""
        )), 1),
        count = c(1, 1, time_length)
      )
    )
  )

  # n(t) is calculated for each monomer
  assign(paste("integral"), cbind(c(1:(time_length)), cumsum(get(
    paste("dn", monomer, sep = "")
  )[, 2] / get(
    paste("longitud_tubo_", monomer, sep = "")
  ))))
  assign(paste("n", monomer, sep = ""), data.frame(matrix(ncol = 2, nrow = 0)))
  assign(paste("n", monomer, sep = ""), rbind(get(paste("n", monomer, sep = "")), integral))

  tmp <- get(paste("n", monomer, sep = ""))
  colnames(tmp) <- c("time", "n")
  assign(paste("n", monomer, sep = ""), tmp)
}

# n(t) is split in in fragments of a defined lenght. To choose this length you have to have in mind that you will need 200 fragments in order to get a decent MSD and not less than 50 points to get a fair regression (the first 10 will be discarded). So, the length and the number of fragments are inversely proportional. If you have a long enough trajectory, this will not be a problem.

rows_n <- nrow(nA)
length_n <-
  50 # length of n(t) fragment. We used this for a 10000 frames trajectory, so we got 200 fragmets.
drop_first_points <- 10

for (monomer in monomer_names) {
  assign(paste("n_frag_", monomer, sep = ""), split(
    get(paste("n", monomer, sep = "")),
    rep(
      1:ceiling(rows_n / length_n),
      each = length_n,
      length.out = rows_n
    )
  ))
  # n2(t) is computed for every fragment and set to n2(0)=0
  n_frag_monomer <- get(paste("n_frag_", monomer, sep = ""))
  a <- 0
  for (i in 1:length(n_frag_monomer)) {
    n_frag_monomer[[i]][, 2] <-
      (n_frag_monomer[[i]][, 2] - n_frag_monomer[[i]][1, 2])^2
    n_frag_monomer[[i]][, 1] <-
      n_frag_monomer[[i]][, 1] - n_frag_monomer[[i]][1, 1]
    a <- n_frag_monomer[[i]][, 2] + a
  }
  # <n2(t)> is computed for every monomer
  b <- cbind(c(1:length(a)), a / (length(n_frag_monomer)))
  assign(paste("MSD_", monomer, sep = ""), data.frame(matrix(ncol = 2, nrow = 0)))
  assign(paste("MSD_", monomer, sep = ""), rbind(get(paste(
    "MSD_", monomer,
    sep = ""
  )), b))

  tmp <- get(paste("MSD_", monomer, sep = ""))
  colnames(tmp) <- c("time", "MSD")
  assign(paste("MSD_", monomer, sep = ""), tmp)

  # the first 10 points of every set are discarded. Refer to Zhu et al. 2004 "It can be shown that when t is much longer than ..."
  assign(paste("MSD_", monomer, "_trunc", sep = ""), data.frame(matrix(ncol = 2, nrow = 0)))
  assign(paste("MSD_", monomer, "_trunc", sep = ""), rbind(get(
    paste("MSD_", monomer, "_trunc", sep = "")
  ), cbind(get(
    paste("MSD_", monomer, sep = "")
  )[drop_first_points:length_n, 1], get(
    paste("MSD_", monomer, sep = "")
  )[drop_first_points:length_n, 2])))

  tmp <- get(paste("MSD_", monomer, "_trunc", sep = ""))
  tmp[, 1] <- tmp[, 1] * frame_step
  assign(paste("MSD_", monomer, "_trunc", sep = ""), tmp)
  # write MSD of each monomer
  assign(
    paste("tname"),
    paste("./", "MSD_", monomer, ".txt", sep = "")
  )

  write.table(
    get(
      paste("MSD_", monomer, "_trunc", sep = "")
    ),
    tname,
    col.names = FALSE,
    row.names = FALSE,
    sep = " "
  )
}

