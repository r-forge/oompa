####################################
# Kelly
kelly <- c(
  white = "#f2f3f4",
  black = "#222222",
  yellow = "#f3c300",
  purple = "#875692",
  orange = "#f38400",
  lightblue = "#a1caf1",
  red = "#be0032",
  buff = "#c2b280",
  gray = "#848482",
  green = "#008856",
  purplishpink = "#e68fac",
  blue = "#0067a5",
  yellowishpink = "#f99379",
  violet = "#604e97",
  orangeyellow = "#f6a600",
  purplishred = "#b3446c",
  greenishyellow = "#dcd300",
  reddishbrown = "#882d17",
  yellowgreen = "#8db600",
  yellowishbrown = "#654522",
  reddishorange = "#e25822",
  olivegreen = "#2b3d26"
  )

####################################
# Green-Armytage

colsch <- list(
  amethyst=c(240, 163, 255),
  blue=c(0, 117, 220),
  caramel=c(153,63,0),
  damson=c(76, 0, 92),
  ebony=c(25, 25, 25),
  forest=c(0, 92, 49),
  green=c(43, 206, 72),
  honeydew=c(255, 204, 153),
  iron=c(128, 128, 128),
  jade=c(148, 255, 181),
  khaki=c(143, 124, 0),
  lime=c(157, 204, 0),
  mallow=c(194, 0, 136),
  navy=c(0, 51, 128),
  orpiment=c(25, 164, 5),
  pink=c(255, 168, 187),
  quagmire=c(66, 102, 0),
  red=c(255, 0, 16),
  sky=c(94, 241, 242),
  turquoise=c(0, 153, 143),
  uranium=c(224, 255, 102),
  violet=c(16, 10, 255),
  wine=c(153, 0, 0),
  xanthin=c(255, 255, 128),
  yellow=c(255, 225, 0),
  zinnia=c(255, 80, 0))
R <- unlist(lapply(colsch, function(x) x[1]/255))
G <- unlist(lapply(colsch, function(x) x[2]/255))
B <- unlist(lapply(colsch, function(x) x[3]/255))
alpha <- as(RGB(R, G, B), "LUV")
green.armytage <- hex(alpha)
names(green.armytage) <- names(colsch)

####################################
# alphacolors

ebony <- hex(LUV(30,0,0))
iron <- hex(LUV(90, 0, 0))
red <- hex(RGB(0.8, 0, 0))
set.seed(567629)
alphacolors <- createPalette(26, c(ebony, iron, red))
d3 <- function(y0) {
  temp <- sweep(alpha@coords, 2, y0, "-")
  dist <- apply(temp^2, 1, sum)
  which(dist==min(dist))
}

alp <- as(hex2RGB(alphacolors), "LUV")
d3(alp@coords[1,])

matcher <- apply(alp@coords, 1, d3)
m <- matcher[1:26]
odd <- which(!(1:26 %in% m))
m[duplicated(m)] <- odd

alphabet <- alphacolors[1:26]
names(alphabet) <- names(colsch)[m]
b <- which(names(alphabet)=='blue')
v <- which(names(alphabet)=='violet')
tu <- which(names(alphabet)=='turquoise')
names(alphabet)[b] <- "turquoise"
names(alphabet)[tu] <- "violet"
names(alphabet)[v] <- "blue"
rm(b, v, tu)
s <- which(names(alphabet) == "sky")
g <- which(names(alphabet) == "green")
j <- which(names(alphabet) == "jade")
k <- which(names(alphabet) == "khaki")
o <- which(names(alphabet) == "orpiment")
w <- which(names(alphabet) == "wine")
names(alphabet)[g] <- "sea"
names(alphabet)[k] <- "jade"
names(alphabet)[o] <- "green"
names(alphabet)[s] <- "wine"
names(alphabet)[j] <- "O"
names(alphabet)[w] <- "kingcrab"
rm(s,g,o,k, j)
names(alphabet)[names(alphabet)=="zinnia"] <- "orange"
names(alphabet)[names(alphabet)=="lime"] <- "zinnia"
names(alphabet)[names(alphabet)=="amethyst"] <- "lavender"
names(alphabet)[names(alphabet)=="O"] <- "amethyst"
names(alphabet)[names(alphabet)=="mallow"] <- "magenta"
names(alphabet)[names(alphabet)=="damson"] <- "ultraviolet"
names(alphabet)[names(alphabet)=="uranium"] <- "damson"
x <- which(names(alphabet) == "xanthin")
h <- which(names(alphabet) == "honeydew")
names(alphabet)[x] <- "honey"
names(alphabet)[h] <- "xanthin"
rm(x,h)
alphabet["ebony"] <- "#565656"
alphabet["iron"] <- "#E2E2E2"
