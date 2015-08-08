{-# LANGUAGE OverloadedStrings #-}

import           Control.Applicative
import           Control.Lens         hiding ((.=))
import           Control.Monad
import           Data.Aeson
import           Data.Aeson.Lens
import qualified Data.ByteString      as S
import qualified Data.ByteString.Lazy as L
import           Text.Printf

readProblem :: Int -> IO Value
readProblem pid = do
    let fname = printf "problems/problem_%d.json" pid
    Just doc <- decode <$> L.readFile fname
    return doc

genAns :: Int -> Int -> Value
genAns pid seed =
    object
        [ "problemId" .= pid
        , "seed"      .= seed
        , "tag"       .= ("´･_･`)(´･_･`" :: String)
        , "solution"  .= (take 100000 $ cycle "Ei!" :: String)
        ]

main :: IO ()
main = do
    rs <- forM [0..23] $ \i -> do
        v <- readProblem i
        let seeds = v^..key"sourceSeeds".values._Integer :: [Integer]
        return $ map (genAns i . fromIntegral) seeds

    L.putStrLn $ encode (concat rs :: [Value])
