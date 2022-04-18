#our base image
FROM julia:1.0.3

# install julia dependencies
COPY deps.jl /usr/src/app/
RUN julia /usr/src/app/deps.jl

# copy files required for the app to run
COPY app.jl /usr/src/app/

# run the application
CMD ["julia", "/usr/src/app/app.jl"]


