library(shiny)
library(lpSolve)

parse_vector <- function(txt, n = 12) {
  x <- as.numeric(trimws(unlist(strsplit(txt, ","))))
  if (length(x) != n || any(is.na(x))) return(NULL)
  x
}

crear_indices <- function(TT) {
  base <- c(P = 0, E = TT, S = 2 * TT, I = 3 * TT, W = 4 * TT, H = 5 * TT, F = 6 * TT)
  function(var, t) base[var] + t
}

extraer_solucion <- function(sol, par) {
  TT <- par$TT
  idx <- crear_indices(TT)

  P <- sol[idx("P", 1:TT)]
  E <- sol[idx("E", 1:TT)]
  S <- sol[idx("S", 1:TT)]
  I <- sol[idx("I", 1:TT)]
  W <- sol[idx("W", 1:TT)]
  H <- sol[idx("H", 1:TT)]
  F <- sol[idx("F", 1:TT)]

  costo_total <- sum(
    par$CN * P + par$CE * E + par$CS * S + par$CI * I +
      par$CW * W + par$CC * H + par$CD * F
  )

  inventario_acumulado <- sum(I)

  detalle <- data.frame(
    Mes = 1:TT,
    Demanda = par$D,
    Produccion_normal = P,
    Produccion_extra = E,
    Subcontratacion = S,
    Inventario_final = I,
    Trabajadores = W,
    Contrataciones = H,
    Despidos = F,
    check.names = FALSE
  )

  list(
    P = P, E = E, S = S, I = I, W = W, H = H, F = F,
    costo_total = costo_total,
    inventario_acumulado = inventario_acumulado,
    detalle = detalle
  )
}

armar_modelo <- function(par, objetivo = c("costo", "inventario"), inv_cap = NULL) {
  objetivo <- match.arg(objetivo)

  TT <- par$TT
  nvars <- 7 * TT
  idx <- crear_indices(TT)

  A <- list()
  dir <- c()
  rhs <- c()

  agregar <- function(row, d, r) {
    A[[length(A) + 1]] <<- row
    dir <<- c(dir, d)
    rhs <<- c(rhs, r)
  }

  for (t in 1:TT) {
    row <- rep(0, nvars)
    row[idx("P", t)] <- 1
    row[idx("E", t)] <- 1
    row[idx("S", t)] <- 1
    row[idx("I", t)] <- -1

    if (t == 1) {
      agregar(row, "=", par$D[t] - par$I0)
    } else {
      row[idx("I", t - 1)] <- 1
      agregar(row, "=", par$D[t])
    }
  }

  for (t in 1:TT) {
    row <- rep(0, nvars)
    row[idx("W", t)] <- 1
    row[idx("H", t)] <- -1
    row[idx("F", t)] <- 1

    if (t == 1) {
      agregar(row, "=", par$W0)
    } else {
      row[idx("W", t - 1)] <- -1
      agregar(row, "=", 0)
    }
  }

  for (t in 1:TT) {
    row <- rep(0, nvars)
    row[idx("P", t)] <- 1
    row[idx("W", t)] <- -par$a
    agregar(row, "<=", 0)
  }

  for (t in 1:TT) {
    row <- rep(0, nvars)
    row[idx("E", t)] <- 1
    agregar(row, "<=", par$HEMAX)
  }

  for (t in 1:TT) {
    row <- rep(0, nvars)
    row[idx("S", t)] <- 1
    agregar(row, "<=", par$SMAX)
  }

  row <- rep(0, nvars)
  row[idx("I", TT)] <- 1
  agregar(row, ">=", par$IFIN)

  if (!is.null(inv_cap)) {
    row <- rep(0, nvars)
    row[idx("I", 1:TT)] <- 1
    agregar(row, "<=", inv_cap)
  }

  mat <- do.call(rbind, A)

  obj <- rep(0, nvars)
  if (objetivo == "costo") {
    obj[idx("P", 1:TT)] <- par$CN
    obj[idx("E", 1:TT)] <- par$CE
    obj[idx("S", 1:TT)] <- par$CS
    obj[idx("I", 1:TT)] <- par$CI
    obj[idx("W", 1:TT)] <- par$CW
    obj[idx("H", 1:TT)] <- par$CC
    obj[idx("F", 1:TT)] <- par$CD
  } else {
    obj[idx("I", 1:TT)] <- 1
  }

  list(mat = mat, dir = dir, rhs = rhs, obj = obj)
}

resolver_modelo <- function(par, objetivo = c("costo", "inventario"), inv_cap = NULL) {
  objetivo <- match.arg(objetivo)
  modelo <- armar_modelo(par, objetivo, inv_cap)

  res <- lp(
    direction = "min",
    objective.in = modelo$obj,
    const.mat = modelo$mat,
    const.dir = modelo$dir,
    const.rhs = modelo$rhs,
    all.int = TRUE
  )

  if (res$status != 0) return(NULL)

  sol <- round(res$solution)
  ext <- extraer_solucion(sol, par)

  list(
    solution = sol,
    obj = res$objval,
    costo_total = ext$costo_total,
    inventario_acumulado = ext$inventario_acumulado,
    detalle = ext$detalle
  )
}

no_dominadas <- function(df) {
  if (nrow(df) <= 1) return(df)

  dom <- rep(FALSE, nrow(df))
  for (i in 1:nrow(df)) {
    for (j in 1:nrow(df)) {
      if (i != j) {
        if (
          df$Costo_total[j] <= df$Costo_total[i] &&
          df$Inventario_acumulado[j] <= df$Inventario_acumulado[i] &&
          (
            df$Costo_total[j] < df$Costo_total[i] ||
            df$Inventario_acumulado[j] < df$Inventario_acumulado[i]
          )
        ) {
          dom[i] <- TRUE
          break
        }
      }
    }
  }
  df[!dom, , drop = FALSE]
}

seleccionar_recomendada <- function(df, criterio = "Equilibrio Pareto") {
  if (nrow(df) == 1) return(df[1, , drop = FALSE])

  if (criterio == "Menor costo") {
    return(df[order(df$Costo_total, df$Inventario_acumulado), , drop = FALSE][1, , drop = FALSE])
  }

  if (criterio == "Menor inventario") {
    return(df[order(df$Inventario_acumulado, df$Costo_total), , drop = FALSE][1, , drop = FALSE])
  }

  cmin <- min(df$Costo_total)
  cmax <- max(df$Costo_total)
  imin <- min(df$Inventario_acumulado)
  imax <- max(df$Inventario_acumulado)

  if (cmax == cmin) {
    df$Costo_norm <- 0
  } else {
    df$Costo_norm <- (df$Costo_total - cmin) / (cmax - cmin)
  }

  if (imax == imin) {
    df$Inv_norm <- 0
  } else {
    df$Inv_norm <- (df$Inventario_acumulado - imin) / (imax - imin)
  }

  df$Distancia_ideal <- sqrt(df$Costo_norm^2 + df$Inv_norm^2)
  df[order(df$Distancia_ideal, df$Costo_total, df$Inventario_acumulado), , drop = FALSE][1, , drop = FALSE]
}

generar_frontera <- function(par, paso = 5) {
  extremo_inv <- resolver_modelo(par, "inventario", NULL)
  extremo_costo <- resolver_modelo(par, "costo", NULL)

  if (is.null(extremo_inv) || is.null(extremo_costo)) return(NULL)

  inv_min <- extremo_inv$inventario_acumulado
  inv_max <- extremo_costo$inventario_acumulado

  caps <- unique(c(seq(inv_min, inv_max, by = max(1, paso)), inv_max))

  filas <- list()
  detalles <- list()
  k <- 1

  for (eps in caps) {
    r <- resolver_modelo(par, "costo", inv_cap = eps)
    if (!is.null(r)) {
      filas[[k]] <- data.frame(
        Solucion = k,
        Limite_inventario = eps,
        Costo_total = r$costo_total,
        Inventario_acumulado = r$inventario_acumulado,
        check.names = FALSE
      )
      detalles[[k]] <- r$detalle
      k <- k + 1
    }
  }

  if (length(filas) == 0) return(NULL)

  frontera <- do.call(rbind, filas)
  frontera <- unique(frontera)
  frontera <- no_dominadas(frontera)
  frontera <- frontera[order(frontera$Inventario_acumulado, frontera$Costo_total), , drop = FALSE]
  frontera$Solucion <- seq_len(nrow(frontera))

  detalle_final <- vector("list", nrow(frontera))
  for (i in 1:nrow(frontera)) {
    eps <- frontera$Limite_inventario[i]
    r <- resolver_modelo(par, "costo", inv_cap = eps)
    detalle_final[[i]] <- r$detalle
  }

  list(
    frontera = frontera,
    inv_min = inv_min,
    inv_max = inv_max,
    detalle = detalle_final
  )
}

ui <- fluidPage(
  titlePanel("Planeacion agregada multiobjetivo para muebles"),

  sidebarLayout(
    sidebarPanel(
      h4("Demandas mensuales"),
      textAreaInput(
        "demanda_txt",
        "Escribe 12 valores separados por coma",
        value = "120,135,150,165,180,195,210,205,190,170,150,130",
        rows = 3
      ),

      tags$hr(),

      h4("Condiciones iniciales"),
      numericInput("I0", "Inventario inicial", value = 30, min = 0, step = 1),
      numericInput("W0", "Trabajadores iniciales", value = 8, min = 0, step = 1),
      numericInput("IFIN", "Inventario final minimo", value = 25, min = 0, step = 1),

      tags$hr(),

      h4("Capacidad"),
      numericInput("a", "Unidades por trabajador por mes", value = 20, min = 1, step = 1),
      numericInput("HEMAX", "Capacidad mensual de tiempo extra", value = 35, min = 0, step = 1),
      numericInput("SMAX", "Capacidad mensual de subcontratacion", value = 25, min = 0, step = 1),

      tags$hr(),

      h4("Costos"),
      numericInput("CN", "Costo produccion normal", value = 1150, min = 0, step = 1),
      numericInput("CE", "Costo produccion extra", value = 1380, min = 0, step = 1),
      numericInput("CS", "Costo subcontratacion", value = 1500, min = 0, step = 1),
      numericInput("CI", "Costo de inventario", value = 45, min = 0, step = 1),
      numericInput("CC", "Costo de contratacion", value = 3500, min = 0, step = 1),
      numericInput("CD", "Costo de despido", value = 4500, min = 0, step = 1),
      numericInput("CW", "Salario mensual por trabajador", value = 9500, min = 0, step = 1),

      tags$hr(),

      numericInput("paso_eps", "Paso para barrer epsilon", value = 5, min = 1, step = 1),
      selectInput(
        "criterio",
        "Criterio para elegir la solucion recomendada",
        choices = c("Equilibrio Pareto", "Menor costo", "Menor inventario"),
        selected = "Equilibrio Pareto"
      ),
      actionButton("resolver", "Calcular frontera de Pareto")
    ),

    mainPanel(
      h3("Resumen"),
      verbatimTextOutput("estado"),

      h3("Datos del modelo"),
      tableOutput("tabla_datos"),

      h3("Frontera de Pareto"),
      tableOutput("tabla_frontera"),

      h3("Plan recomendado"),
      tableOutput("tabla_recomendada"),

      h3("Detalle mensual del plan recomendado"),
      tableOutput("tabla_mensual"),

      h3("Grafica de la frontera"),
      plotOutput("plot_frontera", height = "500px")
    )
  )
)

server <- function(input, output, session) {

  datos_modelo <- eventReactive(input$resolver, {
    D <- parse_vector(input$demanda_txt, 12)

    if (is.null(D)) {
      return(list(error = "La demanda debe tener exactamente 12 valores numericos separados por coma."))
    }

    par <- list(
      TT = 12,
      D = D,
      I0 = input$I0,
      W0 = input$W0,
      IFIN = input$IFIN,
      a = input$a,
      HEMAX = input$HEMAX,
      SMAX = input$SMAX,
      CN = input$CN,
      CE = input$CE,
      CS = input$CS,
      CI = input$CI,
      CC = input$CC,
      CD = input$CD,
      CW = input$CW
    )

    frontera_obj <- generar_frontera(par, input$paso_eps)

    if (is.null(frontera_obj)) {
      return(list(error = "No se encontro solucion factible para esos datos."))
    }

    frontera <- frontera_obj$frontera
    recomendada <- seleccionar_recomendada(frontera, input$criterio)
    idx_rec <- recomendada$Solucion[1]
    detalle_rec <- frontera_obj$detalle[[idx_rec]]

    list(
      error = NULL,
      par = par,
      frontera = frontera,
      recomendada = recomendada,
      detalle_rec = detalle_rec,
      inv_min = frontera_obj$inv_min,
      inv_max = frontera_obj$inv_max
    )
  }, ignoreInit = FALSE)

  output$estado <- renderText({
    d <- datos_modelo()

    if (!is.null(d$error)) return(d$error)

    s <- d$recomendada

    paste0(
      "Soluciones eficientes encontradas: ", nrow(d$frontera), "\n\n",
      "Criterio de seleccion: ", input$criterio, "\n",
      "Solucion recomendada: ", s$Solucion, "\n",
      "Costo total = ", format(round(s$Costo_total, 2), big.mark = ","), "\n",
      "Inventario acumulado = ", format(round(s$Inventario_acumulado, 2), big.mark = ","), "\n",
      "Rango de inventario explorado = [", d$inv_min, ", ", d$inv_max, "]"
    )
  })

  output$tabla_datos <- renderTable({
    d <- datos_modelo()
    validate(need(is.null(d$error), d$error))

    data.frame(
      Concepto = c(
        "Demanda mensual",
        "Inventario inicial",
        "Trabajadores iniciales",
        "Inventario final minimo",
        "Unidades por trabajador",
        "Tiempo extra maximo",
        "Subcontratacion maxima",
        "Costo produccion normal",
        "Costo produccion extra",
        "Costo subcontratacion",
        "Costo inventario",
        "Costo contratacion",
        "Costo despido",
        "Salario mensual"
      ),
      Valor = c(
        paste(d$par$D, collapse = ", "),
        d$par$I0,
        d$par$W0,
        d$par$IFIN,
        d$par$a,
        d$par$HEMAX,
        d$par$SMAX,
        d$par$CN,
        d$par$CE,
        d$par$CS,
        d$par$CI,
        d$par$CC,
        d$par$CD,
        d$par$CW
      ),
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE)

  output$tabla_frontera <- renderTable({
    d <- datos_modelo()
    validate(need(is.null(d$error), d$error))

    df <- d$frontera
    df$Recomendada <- ifelse(df$Solucion == d$recomendada$Solucion[1], "Si", "No")
    df
  }, striped = TRUE, bordered = TRUE, digits = 2)

  output$tabla_recomendada <- renderTable({
    d <- datos_modelo()
    validate(need(is.null(d$error), d$error))

    s <- d$recomendada

    data.frame(
      Indicador = c("Solucion", "Costo total", "Inventario acumulado"),
      Valor = c(
        s$Solucion,
        format(round(s$Costo_total, 2), big.mark = ","),
        format(round(s$Inventario_acumulado, 2), big.mark = ",")
      ),
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE)

  output$tabla_mensual <- renderTable({
    d <- datos_modelo()
    validate(need(is.null(d$error), d$error))
    d$detalle_rec
  }, striped = TRUE, bordered = TRUE, digits = 2)

  output$plot_frontera <- renderPlot({
    d <- datos_modelo()
    validate(need(is.null(d$error), d$error))

    fr <- d$frontera
    s <- d$recomendada

    plot(
      fr$Inventario_acumulado,
      fr$Costo_total,
      type = if (nrow(fr) > 1) "b" else "p",
      pch = 1,
      xlab = "Inventario acumulado",
      ylab = "Costo total",
      main = "Frontera de Pareto"
    )

    text(
      fr$Inventario_acumulado,
      fr$Costo_total,
      labels = fr$Solucion,
      pos = 3,
      cex = 0.8
    )

    points(
      s$Inventario_acumulado,
      s$Costo_total,
      pch = 19,
      cex = 1.4
    )
  })
}

shinyApp(ui = ui, server = server)