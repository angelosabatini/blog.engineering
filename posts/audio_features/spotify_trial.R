library(dplyr)
library(spotifyr)
library(plotly)
library(ggplot2)

id     <- "fd6b24ef562746b3a4fc04776d9581f4"
secret <- "02af4b54943b4ccc806a0ce5127917fa"

Sys.setenv(SPOTIFY_CLIENT_ID = id)
Sys.setenv(SPOTIFY_CLIENT_SECRET = secret)

access_token <- get_spotify_access_token()

my_id     <- "lzdooepemuu0wu7n18rlzq40n"
my_plists <- get_user_playlists(my_id)

playlist_id <- my_plists %>%
  filter(name == 'Hit parade italiane (1960-2019)') %>%
  pull(id)

tracks  <- get_playlist_tracks(playlist_id = playlist_id, limit = 100, authorization = access_token)
stringa <- tracks$track.id

popularity <- tracks$track.popularity

features <- get_track_audio_features(stringa)

l_stringa <- length(stringa)
artist <- array(0, l_stringa)
track <- array(0, l_stringa)
popularity <- array(0, l_stringa)
release <- array(0, l_stringa)
album <- array(0, l_stringa)

for (i in 1:l_stringa) {
  artist[i] <- get_track(stringa[i])$artists$name
  track[i]  <- get_track(stringa[i])$name
  popularity[i] <- get_track(stringa[i])$popularity
  release[i] <- get_track(stringa[i])$album$release_date
  album[i] <- get_track(stringa[i])$album$name
}

features <- features %>%
  mutate(track = track, album = album, release = release, popularity = popularity) %>%
  filter(album == "Un Sabato Italiano") %>%
  select(track, popularity, danceability:tempo, duration_ms:time_signature)
